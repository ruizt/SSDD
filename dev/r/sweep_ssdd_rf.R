## sweep_ssdd_rf.R
##
## Radius sensitivity sweep for the 4-metric SSDD RF baseline.
## Iterates over the (r_D, r_S) grid in _data/processed/sweep/sweep_all.csv
## using the same spatial block CV as baseline_ssdd_rf.R.
##
## Fold assignments are computed ONCE per fire (coordinates are constant
## across radii), then reused for every (r_D, r_S) combination.
##
## Output: dev/r/sweep_results.csv — per-fold AUC/log-loss/Brier for every
##         (r_D, r_S, setting, fold) combination.

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

# ── 1. Load and prepare data ─────────────────────────────────────────────────

message("Loading sweep data...")
sweep <- read_csv(
  file.path(PROCESSED, "sweep", "sweep_all.csv"),
  show_col_types = FALSE
)

# Nearest-neighbor proximity join with kenny covariates (within 5 m).
# Done once per fire on a reference radius combo; matched ssdd_id set is
# then applied to all (r_D, r_S) combos.
join_fire <- function(fire, sweep_data) {
  kenny <- read_csv(
    file.path(PROCESSED, fire, "covariates",
              sprintf("%s_burned_struc_model_inputs.csv", fire)),
    show_col_types = FALSE
  )

  # Use one (r_D, r_S) slice for the proximity join — coords are constant
  first_combo <- sweep_data |> distinct(r_D, r_S) |> slice(1)

  ref <- sweep_data |>
    filter(fire == !!fire,
           r_D == first_combo$r_D, r_S == first_combo$r_S,
           !is.na(DAMAGE))

  ssdd_sf  <- st_as_sf(ref,   coords = c("cent_x", "cent_y"), crs = EPSG)
  kenny_sf <- st_as_sf(kenny, coords = c("utm_x",  "utm_y"),  crs = EPSG)

  nn_idx  <- st_nearest_feature(ssdd_sf, kenny_sf)
  nn_dist <- as.numeric(st_distance(ssdd_sf, kenny_sf[nn_idx, ], by_element = TRUE))
  matched <- nn_dist <= 5

  matched_ids <- ref$ssdd_id[matched]
  terrain_df  <- kenny[nn_idx[matched], c("UID", TERRAIN_COLS)]

  # Map ssdd_id -> terrain for this fire
  id_terrain <- bind_cols(tibble(ssdd_id = matched_ids), terrain_df)

  # Apply to all (r_D, r_S) combos
  sweep_data |>
    filter(fire == !!fire, !is.na(DAMAGE), ssdd_id %in% matched_ids) |>
    mutate(destroyed = as.integer(DAMAGE == "Destroyed (>50%)")) |>
    left_join(id_terrain, by = "ssdd_id") |>
    select(fire, r_D, r_S, cent_x, cent_y,
           all_of(SSDD_COLS), all_of(TERRAIN_COLS), destroyed)
}

sweep_joined <- bind_rows(lapply(FIRES, join_fire, sweep_data = sweep))

message(sprintf("Sweep joined: %d rows", nrow(sweep_joined)))

# Check sample sizes per fire (should be constant across radii)
sweep_joined |>
  filter(r_D == min(r_D), r_S == min(r_S)) |>
  count(fire) |>
  print()

# ── 2. Compute fold assignments once per fire ─────────────────────────────────

# Coordinates are the same for every (r_D, r_S), so fold IDs are reusable.
first_combo <- sweep_joined |> distinct(r_D, r_S) |> slice(1)
ref <- sweep_joined |> filter(r_D == first_combo$r_D, r_S == first_combo$r_S)

set.seed(7291)
fold_ids_list <- list()
for (f in FIRES) {
  d <- ref |> filter(fire == f)
  sf_d <- st_as_sf(d, coords = c("cent_x", "cent_y"), crs = EPSG, remove = FALSE)
  fold_obj <- cv_spatial(
    x = sf_d, k = N_FOLDS, size = BLOCK_SIZE,
    selection = "random", iteration = 50, seed = 7291,
    progress = FALSE, plot = FALSE, report = FALSE
  )
  fold_ids_list[[f]] <- fold_obj$folds_ids
  message(sprintf("[%s] %d buildings, fold sizes: %s",
                  f, nrow(d),
                  paste(table(fold_obj$folds_ids), collapse = "/")))
}

# ── 3. Scoring helpers ───────────────────────────────────────────────────────

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
  tb <- tibble(
    truth = factor(truth, levels = c(0L, 1L)),
    pred  = pred
  )
  tibble(
    auc      = roc_auc_vec(tb$truth, tb$pred, event_level = "second"),
    log_loss = mn_log_loss_vec(tb$truth, tb$pred, event_level = "second"),
    brier    = brier_class_vec(tb$truth, tb$pred, event_level = "second")
  )
}

# ── 4. Run sweep ─────────────────────────────────────────────────────────────

radius_grid <- sweep_joined |> distinct(r_D, r_S) |> arrange(r_D, r_S)
message(sprintf("Running %d radius combinations...", nrow(radius_grid)))

all_sweep_results <- list()
counter <- 0L

for (i in seq_len(nrow(radius_grid))) {
  rD <- radius_grid$r_D[i]
  rS <- radius_grid$r_S[i]

  sub <- sweep_joined |> filter(r_D == rD, r_S == rS)
  per_fire_sub <- setNames(
    lapply(FIRES, \(f) sub |> filter(fire == f)),
    FIRES
  )
  pooled_sub <- sub

  fold_results <- list()

  # Within-fire block CV for each fire
  for (f in FIRES) {
    d <- per_fire_sub[[f]]
    fids <- fold_ids_list[[f]]

    for (k in seq_len(N_FOLDS)) {
      out <- fit_rf(d[fids != k, ], d[fids == k, ])
      sc  <- score_fold(out$truth, out$pred)
      fold_results[[length(fold_results) + 1]] <- tibble(
        r_D = rD, r_S = rS,
        setting = paste0(f, "-within"),
        fold = as.character(k),
        auc = sc$auc, log_loss = sc$log_loss, brier = sc$brier
      )
    }
  }

  # Pooled within-fire block CV
  pooled_fids <- c(fold_ids_list$eaton, fold_ids_list$palisades)
  for (k in seq_len(N_FOLDS)) {
    out <- fit_rf(pooled_sub[pooled_fids != k, ], pooled_sub[pooled_fids == k, ])
    sc  <- score_fold(out$truth, out$pred)
    fold_results[[length(fold_results) + 1]] <- tibble(
      r_D = rD, r_S = rS,
      setting = "pooled-within", fold = as.character(k),
      auc = sc$auc, log_loss = sc$log_loss, brier = sc$brier
    )
  }

  # LOFO
  for (test_fire in FIRES) {
    out <- fit_rf(
      pooled_sub |> filter(fire != test_fire),
      pooled_sub |> filter(fire == test_fire)
    )
    sc <- score_fold(out$truth, out$pred)
    fold_results[[length(fold_results) + 1]] <- tibble(
      r_D = rD, r_S = rS,
      setting = "lofo", fold = test_fire,
      auc = sc$auc, log_loss = sc$log_loss, brier = sc$brier
    )
  }

  all_sweep_results[[i]] <- bind_rows(fold_results)
  counter <- counter + 1L
  message(sprintf("  [%d/%d] rD=%d rS=%d done", counter, nrow(radius_grid), rD, rS))
}

sweep_results <- bind_rows(all_sweep_results)

out_path <- "dev/r/sweep_results.csv"
write_csv(sweep_results, out_path)
message(sprintf("Wrote %s (%d rows)", out_path, nrow(sweep_results)))

# ── 5. Summary ───────────────────────────────────────────────────────────────

sweep_summary <- sweep_results |>
  group_by(r_D, r_S, setting) |>
  summarise(
    auc_mean = mean(auc, na.rm = TRUE),
    auc_sd   = sd(auc, na.rm = TRUE),
    .groups  = "drop"
  )

sweep_summary |> print(n = Inf)
