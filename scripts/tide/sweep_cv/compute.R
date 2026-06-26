#!/usr/bin/env Rscript
## compute.R — per-(r_D, r_S) CV entrypoint for the Tide CV sweep.
##
## Reads sweep_all.csv from the data PVC, filters to ONE (r_D, r_S) combo,
## NN-joins Kenny terrain covariates, computes spatial-block + LOFO folds,
## and runs all CV fits in parallel via parallel::mclapply (one fork per
## fit, ranger single-threaded inside each fork).
##
## Env vars
## --------
## SSDD_R_D        SD radius (m) — required
## SSDD_R_S        SS radius (m) — required
## SSDD_CORES      Forks for mclapply (default 4)
## SSDD_N_FOLDS    Within-fire block CV folds (default 10)
## SSDD_BLOCK_SIZE Spatial block size in meters (default 500)
## SSDD_DATA_DIR   Input root (default /data); expects:
##                   <data>/sweep/sweep_all.csv
##                   <data>/<fire>/covariates/<fire>_burned_struc_model_inputs.csv
## SSDD_OUT_DIR    Output root (default /jobs/output)
##
## Output
## ------
##   <out>/rD<r_D>_rS<r_S>/rD<r_D>_rS<r_S>_cv.csv
##   columns: r_D, r_S, setting, fold, auc, log_loss, brier

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(blockCV)
  library(ranger)
  library(yardstick)
  library(parallel)
})

# ----- Env --------------------------------------------------------------------

rD          <- as.integer(Sys.getenv("SSDD_R_D"))
rS          <- as.integer(Sys.getenv("SSDD_R_S"))
n_cores     <- as.integer(Sys.getenv("SSDD_CORES", "4"))
N_FOLDS     <- as.integer(Sys.getenv("SSDD_N_FOLDS", "10"))
BLOCK_SIZE  <- as.integer(Sys.getenv("SSDD_BLOCK_SIZE", "500"))
data_dir    <- Sys.getenv("SSDD_DATA_DIR", "/data")
out_dir     <- Sys.getenv("SSDD_OUT_DIR",  "/jobs/output")

stopifnot(!is.na(rD), rD > 0L, !is.na(rS), rS > 0L, n_cores > 0L, N_FOLDS > 1L)

FIRES        <- c("eaton", "palisades")
EPSG         <- 32611L
SSDD_COLS    <- c("KD_raw", "BA_raw", "DP_raw", "OP_raw")
TERRAIN_COLS <- c("elevation", "slope", "aspect")
FEATURE_COLS <- c(SSDD_COLS, TERRAIN_COLS)

message(sprintf("[cv-sweep] rD=%d  rS=%d  cores=%d  folds=%d  block=%dm",
                rD, rS, n_cores, N_FOLDS, BLOCK_SIZE))

# ----- Data load + per-fire NN join with Kenny terrain ------------------------

sweep_path <- file.path(data_dir, "sweep", "sweep_all.csv")
message(sprintf("[cv-sweep] reading %s", sweep_path))
sweep <- read_csv(sweep_path, show_col_types = FALSE) |>
  filter(r_D == rD, r_S == rS)
message(sprintf("[cv-sweep] %d rows for this combo before NN join", nrow(sweep)))

join_fire <- function(fire) {
  kenny_path <- file.path(
    data_dir, fire, "covariates",
    sprintf("%s_burned_struc_model_inputs.csv", fire)
  )
  kenny <- read_csv(kenny_path, show_col_types = FALSE)

  ref <- sweep |> filter(fire == !!fire, !is.na(DAMAGE))
  ssdd_sf  <- st_as_sf(ref,   coords = c("cent_x", "cent_y"), crs = EPSG)
  kenny_sf <- st_as_sf(kenny, coords = c("utm_x",  "utm_y"),  crs = EPSG)

  nn_idx  <- st_nearest_feature(ssdd_sf, kenny_sf)
  nn_dist <- as.numeric(st_distance(ssdd_sf, kenny_sf[nn_idx, ], by_element = TRUE))
  matched <- nn_dist <= 5

  matched_ids <- ref$ssdd_id[matched]
  id_terrain  <- bind_cols(
    tibble(ssdd_id = matched_ids),
    kenny[nn_idx[matched], TERRAIN_COLS]
  )

  ref |>
    filter(ssdd_id %in% matched_ids) |>
    mutate(destroyed = as.integer(DAMAGE == "Destroyed (>50%)")) |>
    left_join(id_terrain, by = "ssdd_id") |>
    select(fire, cent_x, cent_y, all_of(FEATURE_COLS), destroyed)
}

joined <- bind_rows(lapply(FIRES, join_fire))
message(sprintf("[cv-sweep] %d rows after Kenny join (%s)",
                nrow(joined),
                paste(sprintf("%s=%d", FIRES,
                              vapply(FIRES,
                                     \(f) sum(joined$fire == f),
                                     integer(1))),
                      collapse = ", ")))

# ----- Spatial block CV folds (deterministic per fire) ------------------------

set.seed(7291)
fold_ids_list <- setNames(vector("list", length(FIRES)), FIRES)
for (f in FIRES) {
  d <- joined |> filter(fire == f)
  sf_d <- st_as_sf(d, coords = c("cent_x", "cent_y"), crs = EPSG, remove = FALSE)
  fold_obj <- cv_spatial(
    x         = sf_d,
    k         = N_FOLDS,
    size      = BLOCK_SIZE,
    selection = "random",
    iteration = 50,
    seed      = 7291,
    progress  = FALSE,
    plot      = FALSE,
    report    = FALSE
  )
  fold_ids_list[[f]] <- fold_obj$folds_ids
}

# ----- Build the per-fit work list -------------------------------------------

per_fire_data <- setNames(lapply(FIRES, \(f) joined |> filter(fire == f)), FIRES)

specs <- list()

# Within-fire spatial-block CV
for (f in FIRES) {
  d    <- per_fire_data[[f]]
  fids <- fold_ids_list[[f]]
  for (k in seq_len(N_FOLDS)) {
    specs[[length(specs) + 1]] <- list(
      train   = d[fids != k, ],
      test    = d[fids == k, ],
      setting = paste0(f, "-within"),
      fold    = as.character(k)
    )
  }
}

# Pooled within-fire (fold id is the per-fire fold concatenated by fire order)
pooled_fids <- c(fold_ids_list$eaton, fold_ids_list$palisades)
for (k in seq_len(N_FOLDS)) {
  specs[[length(specs) + 1]] <- list(
    train   = joined[pooled_fids != k, ],
    test    = joined[pooled_fids == k, ],
    setting = "pooled-within",
    fold    = as.character(k)
  )
}

# LOFO (one fit per held-out fire)
for (test_fire in FIRES) {
  specs[[length(specs) + 1]] <- list(
    train   = joined |> filter(fire != test_fire),
    test    = joined |> filter(fire == test_fire),
    setting = "lofo",
    fold    = test_fire
  )
}

message(sprintf("[cv-sweep] %d CV fits queued", length(specs)))

# ----- Per-fit RF + scoring --------------------------------------------------

fit_rf <- function(train, test) {
  needed <- c("destroyed", FEATURE_COLS)
  train  <- train[complete.cases(train[, needed]), , drop = FALSE]
  test   <- test [complete.cases(test [, needed]), , drop = FALSE]

  rf_train <- data.frame(
    destroyed = factor(train$destroyed, levels = c(0L, 1L)),
    train[, FEATURE_COLS, drop = FALSE]
  )
  rf_fit <- ranger(
    destroyed ~ ., data = rf_train,
    probability = TRUE, num.trees = 500, seed = 7291,
    num.threads = 1L  # mclapply provides the outer parallelism
  )
  pred <- predict(rf_fit, data = test[, FEATURE_COLS, drop = FALSE])$predictions[, "1"]
  list(truth = test$destroyed, pred = pred)
}

score_fold <- function(truth, pred) {
  tb <- tibble(truth = factor(truth, levels = c(0L, 1L)), pred = pred)
  tibble(
    auc      = roc_auc_vec     (tb$truth, tb$pred, event_level = "second"),
    log_loss = mn_log_loss_vec (tb$truth, tb$pred, event_level = "second"),
    brier    = brier_class_vec (tb$truth, tb$pred, event_level = "second")
  )
}

# ----- Run fits in parallel --------------------------------------------------

t0 <- Sys.time()
message(sprintf("[cv-sweep] fitting on %d cores ...", n_cores))
results <- mclapply(specs, function(spec) {
  out <- fit_rf(spec$train, spec$test)
  sc  <- score_fold(out$truth, out$pred)
  tibble(
    r_D     = rD,
    r_S     = rS,
    setting = spec$setting,
    fold    = spec$fold,
    auc     = sc$auc,
    log_loss = sc$log_loss,
    brier   = sc$brier
  )
}, mc.cores = n_cores)

# mc.preschedule defaults can swallow errors silently; catch any non-tibble.
errs <- vapply(results, \(x) !inherits(x, "data.frame"), logical(1))
if (any(errs)) {
  stop(sprintf("%d fits failed; check container logs", sum(errs)))
}

results <- bind_rows(results)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
message(sprintf("[cv-sweep] %d fits in %.1f sec (%.2f sec/fit avg)",
                nrow(results), elapsed, elapsed / nrow(results)))

# ----- Write per-job CSV ------------------------------------------------------

run_dir  <- file.path(out_dir, sprintf("rD%d_rS%d", rD, rS))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(run_dir, sprintf("rD%d_rS%d_cv.csv", rD, rS))
write_csv(results, out_path)
message(sprintf("[cv-sweep] wrote %s", out_path))
