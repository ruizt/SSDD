## baseline_ssdd_rf.R
##
## Baseline 1: All four raw SSDD metrics + terrain fed directly into RF.
## Tests whether OP_raw (orientation proximity) adds signal beyond the
## three benchmark-equivalent metrics (KD, BA, DP). No blending — RF
## learns how to use all four metrics freely.
##
## Comparing this to benchmark 2 (3-metric) isolates OP_raw's marginal
## contribution. Comparing to benchmark 1 (build_dens) shows whether the
## full SSDD feature set outperforms simple building density.
##
## Features: KD_raw, BA_raw, DP_raw, OP_raw, elevation, slope, aspect
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
FEATURE_COLS <- c(SSDD_COLS, TERRAIN_COLS)

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

fit_rf <- function(train, test, features = FEATURE_COLS) {
  needed <- c("destroyed", features)
  train  <- train[complete.cases(train[, needed]), , drop = FALSE]
  test   <- test [complete.cases(test [, needed]), , drop = FALSE]

  rf_train <- data.frame(
    destroyed = factor(train$destroyed, levels = c(0L, 1L)),
    train[, features, drop = FALSE]
  )
  rf_fit  <- ranger(destroyed ~ ., data = rf_train,
                    probability = TRUE, num.trees = 500, seed = 7291)
  rf_pred <- predict(rf_fit, data = test[, features, drop = FALSE])$predictions[, "1"]

  list(truth = test$destroyed, pred = rf_pred)
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
  tibble(fold = fold_label, auc = scores$auc,
         log_loss = scores$log_loss, brier = scores$brier)
}

make_fold_preds <- function(fold_label, out) {
  tibble(fold = fold_label, truth = out$truth, pred = out$pred)
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
    out <- fit_rf(data[fold_ids != f, ], data[fold_ids == f, ])
    list(
      scores = make_fold_row(as.character(f), out),
      preds  = make_fold_preds(as.character(f), out)
    )
  })
  list(
    scores = bind_rows(lapply(per_fold, `[[`, "scores")) |> mutate(setting = label, .before = 1),
    preds  = bind_rows(lapply(per_fold, `[[`, "preds"))  |> mutate(setting = label, .before = 1)
  )
}

run_lofo <- function(pooled) {
  message("[pooled] leave-one-fire-out CV")
  fires <- unique(pooled$fire)
  per_fire <- lapply(fires, function(test_fire) {
    out <- fit_rf(
      pooled |> filter(fire != test_fire),
      pooled |> filter(fire == test_fire)
    )
    list(
      scores = make_fold_row(test_fire, out),
      preds  = make_fold_preds(test_fire, out)
    )
  })
  list(
    scores = bind_rows(lapply(per_fire, `[[`, "scores")) |> mutate(setting = "lofo", .before = 1),
    preds  = bind_rows(lapply(per_fire, `[[`, "preds"))  |> mutate(setting = "lofo", .before = 1)
  )
}

summarise_results <- function(results) {
  results |>
    group_by(setting) |>
    summarise(
      auc_mean = mean(auc, na.rm = TRUE),
      auc_sd   = sd(auc, na.rm = TRUE),
      log_loss_mean = mean(log_loss, na.rm = TRUE),
      brier_mean = mean(brier, na.rm = TRUE),
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

res_eaton     <- run_within_fire_cv(per_fire_data$eaton,     "eaton-within")
res_palisades <- run_within_fire_cv(per_fire_data$palisades, "palisades-within")
res_pooled    <- run_within_fire_cv(pooled,                  "pooled-within")
res_lofo      <- run_lofo(pooled)

# ── 7. Combine + summarise ───────────────────────────────────────────────────

results <- bind_rows(res_eaton$scores, res_palisades$scores,
                     res_pooled$scores, res_lofo$scores)
predictions <- bind_rows(res_eaton$preds, res_palisades$preds,
                         res_pooled$preds, res_lofo$preds)

out_path <- "dev/r/baseline_ssdd_rf_results.csv"
write_csv(results, out_path)
message(sprintf("Wrote %s", out_path))

pred_path <- "dev/r/baseline_ssdd_rf_predictions.csv"
write_csv(predictions, pred_path)
message(sprintf("Wrote %s", pred_path))

summarise_results(results)
