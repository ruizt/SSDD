## baseline_blend_rf.R
##
## Baseline 2: Two-step logistic-blend + RF. Produces a single univariate
## SSDD score by learning blend weights over the four raw metrics, then
## feeds that score + terrain into ranger.
##
## For each CV fold:
##   1. Fit logistic regression (destroyed ~ KD + BA + DP + OP) on TRAIN
##      to learn blend weights.
##   2. Compute the linear predictor (sans intercept) as a univariate
##      SSDD blend score for TRAIN and TEST.
##   3. Fit ranger on (blend_score + terrain) on TRAIN; predict TEST.
##   4. Score: AUC, log-loss, Brier.
##
## Purpose: test whether a supervised univariate blend can match the
## performance of using all four metrics separately (baseline 1). The
## blend is the target deliverable — a single SSDD score to replace
## build_dens in downstream models.
##
## Features: blend_score (from KD, BA, DP, OP), elevation, slope, aspect
##
## Inputs:
##   _data/processed/<fire>/<fire>_raw_metrics.csv
##   _data/processed/<fire>/covariates/<fire>_burned_struc_model_inputs.csv  (terrain only)

# ── 0. Libraries + config ────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(blockCV)
  library(ranger)
  library(yardstick)
})

PROCESSED  <- "_data/processed"
FIRES      <- c("eaton", "palisades")
EPSG       <- 32611L
N_FOLDS    <- 10L
BLOCK_SIZE <- 500

SSDD_COLS    <- c("KD_raw", "BA_raw", "DP_raw", "OP_raw")
TERRAIN_COLS <- c("elevation", "slope", "aspect")

# ── 1. Function definitions ──────────────────────────────────────────────────

read_fire <- function(fire) {
  ssdd <- read_csv(
    file.path(PROCESSED, fire, sprintf("%s_raw_metrics.csv", fire)),
    show_col_types = FALSE
  )

  kenny <- read_csv(
    file.path(PROCESSED, fire, "covariates",
              sprintf("%s_burned_struc_model_inputs.csv", fire)),
    show_col_types = FALSE
  )

  # Nearest-neighbor proximity join within 5 m
  ssdd_sf  <- st_as_sf(ssdd,  coords = c("cent_x", "cent_y"), crs = EPSG)
  kenny_sf <- st_as_sf(kenny, coords = c("utm_x",  "utm_y"),  crs = EPSG)

  nn_idx  <- st_nearest_feature(ssdd_sf, kenny_sf)
  nn_dist <- as.numeric(st_distance(ssdd_sf, kenny_sf[nn_idx, ], by_element = TRUE))
  matched <- nn_dist <= 5

  ssdd_keep  <- c("ssdd_id", "cent_x", "cent_y", SSDD_COLS)
  kenny_keep <- c("UID", "DAMAGE", TERRAIN_COLS)

  joined <- bind_cols(
    ssdd[matched, ssdd_keep],
    kenny[nn_idx[matched], kenny_keep]
  )

  joined |>
    mutate(
      fire      = fire,
      destroyed = as.integer(DAMAGE == "Destroyed (>50%)")
    )
}

fit_blend_rf <- function(train, test,
                         ssdd_cols    = SSDD_COLS,
                         terrain_cols = TERRAIN_COLS) {
  needed <- c("destroyed", ssdd_cols, terrain_cols)
  train  <- train[complete.cases(train[, needed]), , drop = FALSE]
  test   <- test [complete.cases(test [, needed]), , drop = FALSE]

  # Step 1: logistic blend — learn combination weights on train
  blend_f   <- as.formula(paste("destroyed ~", paste(ssdd_cols, collapse = " + ")))
  blend_glm <- glm(blend_f, data = train, family = binomial())
  blend_coefs <- coef(blend_glm)[ssdd_cols]

  # Linear predictor (without intercept) as the univariate SSDD score
  blend_train <- as.numeric(as.matrix(train[, ssdd_cols]) %*% blend_coefs)
  blend_test  <- as.numeric(as.matrix(test[, ssdd_cols])  %*% blend_coefs)

  # Step 2: RF on (blend_score + terrain)
  rf_train <- data.frame(
    destroyed   = factor(train$destroyed, levels = c(0L, 1L)),
    blend_score = blend_train,
    train[, terrain_cols, drop = FALSE]
  )
  rf_test <- data.frame(
    blend_score = blend_test,
    test[, terrain_cols, drop = FALSE]
  )

  rf_fit  <- ranger(destroyed ~ ., data = rf_train,
                    probability = TRUE, num.trees = 500, seed = 7291)
  rf_pred <- predict(rf_fit, data = rf_test)$predictions[, "1"]

  list(blend_coefs = blend_coefs,
       truth       = test$destroyed,
       pred        = rf_pred)
}

score_fold <- function(truth, pred) {
  tibble(
    truth = factor(truth, levels = c(0L, 1L)),
    pred  = pred
  ) |>
    summarise(
      auc      = roc_auc_vec(truth, pred, event_level = "second"),
      log_loss = mn_log_loss_vec(truth, pred, event_level = "second"),
      brier    = brier_class_vec(truth, pred, event_level = "second")
    )
}

make_fold_row <- function(fold_label, out) {
  scores <- score_fold(out$truth, out$pred)
  tibble(
    fold = fold_label, auc = scores$auc,
    log_loss = scores$log_loss, brier = scores$brier,
    blend_KD = out$blend_coefs["KD_raw"],
    blend_BA = out$blend_coefs["BA_raw"],
    blend_DP = out$blend_coefs["DP_raw"],
    blend_OP = out$blend_coefs["OP_raw"]
  )
}

run_within_fire_cv <- function(data, label, n_folds = N_FOLDS, block_m = BLOCK_SIZE) {
  message(sprintf("[%s] block CV (k=%d, ~%dm)", label, n_folds, block_m))
  sf_data <- st_as_sf(data, coords = c("cent_x", "cent_y"), crs = EPSG, remove = FALSE)
  fold_obj <- cv_spatial(
    x = sf_data, k = n_folds, size = block_m,
    selection = "random", iteration = 50, seed = 7291,
    progress = FALSE, plot = FALSE, report = FALSE
  )
  fold_ids <- fold_obj$folds_ids

  per_fold <- lapply(seq_len(n_folds), function(f) {
    out <- fit_blend_rf(data[fold_ids != f, ], data[fold_ids == f, ])
    make_fold_row(as.character(f), out)
  })
  bind_rows(per_fold) |> mutate(setting = label, .before = 1)
}

run_lofo <- function(pooled) {
  message("[pooled] leave-one-fire-out CV")
  fires <- unique(pooled$fire)
  per_fire <- lapply(fires, function(test_fire) {
    out <- fit_blend_rf(
      pooled |> filter(fire != test_fire),
      pooled |> filter(fire == test_fire)
    )
    make_fold_row(test_fire, out)
  })
  bind_rows(per_fire) |> mutate(setting = "lofo", .before = 1)
}

summarise_results <- function(results) {
  results |>
    group_by(setting) |>
    summarise(
      auc_mean = mean(auc, na.rm = TRUE),
      auc_sd   = sd(auc, na.rm = TRUE),
      log_loss_mean = mean(log_loss, na.rm = TRUE),
      brier_mean = mean(brier, na.rm = TRUE),
      blend_KD_mean = mean(blend_KD, na.rm = TRUE),
      blend_BA_mean = mean(blend_BA, na.rm = TRUE),
      blend_DP_mean = mean(blend_DP, na.rm = TRUE),
      blend_OP_mean = mean(blend_OP, na.rm = TRUE),
      .groups = "drop"
    )
}

# ── 2. Load data ─────────────────────────────────────────────────────────────

set.seed(7291)
per_fire_data <- setNames(lapply(FIRES, read_fire), FIRES)
for (fire in FIRES) {
  message(sprintf("[%s] N=%d; %d destroyed",
                  fire, nrow(per_fire_data[[fire]]),
                  sum(per_fire_data[[fire]]$destroyed)))
}
pooled <- bind_rows(per_fire_data)

# ── 3. Within-fire block CV: Eaton ────────────────────────────────────────────

res_eaton <- run_within_fire_cv(per_fire_data$eaton, "eaton-within")
res_eaton

# ── 4. Within-fire block CV: Palisades ───────────────────────────────────────

res_palisades <- run_within_fire_cv(per_fire_data$palisades, "palisades-within")
res_palisades

# ── 5. Within-fire block CV: Pooled ──────────────────────────────────────────

res_pooled <- run_within_fire_cv(pooled, "pooled-within")
res_pooled

# ── 6. Leave-one-fire-out CV ─────────────────────────────────────────────────

res_lofo <- run_lofo(pooled)
res_lofo

# ── 7. Combine + summarise ───────────────────────────────────────────────────

results <- bind_rows(res_eaton, res_palisades, res_pooled, res_lofo)

out_path <- "dev/r/baseline_blend_rf_results.csv"
write_csv(results, out_path)
message(sprintf("Wrote %s", out_path))

summarise_results(results)
