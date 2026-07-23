#!/usr/bin/env Rscript
## compute.R — per-r_D CV entrypoint for the Tide CV sweep.
##
## Reads sweep_all.csv (from sweep_process), filters to ONE r_D, NN-joins terrain
## covariates, and evaluates every candidate metric FORM-SETTING at two scales:
##
##   structure scale  — RF classifier: {SD, SS} + terrain -> destroyed (0/1)
##   neighbourhood    — aggregate-then-fit: fine-block means of {SD, SS, terrain}
##                      -> block damage-rate (RF regression)
##
## FORM-SETTINGS (12) = 4 SD columns x 3 SS columns, one {SD, SS} pair per fit.
##
## CV scheme: ONE pooled spatial-block CV (no LOFO). Each fold holds out a
## spatial block from EVERY fire (per-fire blocks, shared fold index), so held-out
## performance can be disaggregated by fire. Fold blocks are COARSE (>= a few x
## the largest r_D); neighbourhood aggregation uses a separate FINE block grid,
## each fine block inheriting its structures' coarse fold.
##
## Secondary: per-fire within-fire CV for the shipped form only, as a
## fire-heterogeneity diagnostic (does a fire predict itself?).
##
## Env vars
## --------
## SSDD_R_D         SD radius (m) — required
## SSDD_CORES       forks for mclapply (default 4)
## SSDD_N_FOLDS     spatial CV folds (default 8)
## SSDD_CV_BLOCK    coarse fold block size, m (default 1500)
## SSDD_NBHD_BLOCK  fine neighbourhood block size, m (default 400)
## SSDD_DATA_DIR    input root (default /data): <data>/sweep/sweep_all.csv,
##                  <data>/<fire>/covariates/<fire>_burned_struc_model_inputs.csv
## SSDD_OUT_DIR     output root (default /jobs/output)
##
## Output: <out>/rD<r_D>/rD<r_D>_cv.csv with columns
##   r_D, scale, cv, form, sd_form, ss_form, fire, fold, n, score, aux
##   scale ∈ {structure, neighbourhood}; cv ∈ {pooled, within}; fire includes
##   "pooled"; score = AUC (structure) / concordance (neighbourhood);
##   aux = log_loss (structure) / rmse (neighbourhood).

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(sf); library(blockCV)
  library(ranger); library(yardstick); library(parallel)
})

# ----- Env --------------------------------------------------------------------
rD         <- as.integer(Sys.getenv("SSDD_R_D"))
n_cores    <- as.integer(Sys.getenv("SSDD_CORES", "4"))
N_FOLDS    <- as.integer(Sys.getenv("SSDD_N_FOLDS", "8"))
CV_BLOCK   <- as.integer(Sys.getenv("SSDD_CV_BLOCK", "1500"))
NBHD_BLOCK <- as.integer(Sys.getenv("SSDD_NBHD_BLOCK", "400"))
data_dir   <- Sys.getenv("SSDD_DATA_DIR", "/data")
out_dir    <- Sys.getenv("SSDD_OUT_DIR",  "/jobs/output")
stopifnot(!is.na(rD), rD > 0L, n_cores > 0L, N_FOLDS > 1L)

FIRES        <- c("eaton", "palisades", "mountain")
EPSG         <- 32611L
SD_FORMS     <- c("SD_uniform_root_area", "SD_uniform_unit",
                  "SD_quartic_root_area", "SD_quartic_unit")
SS_FORMS     <- c("SS_flat", "SS_gauss", "SS_cos2")
TERRAIN      <- c("elevation", "slope", "aspect")
SHIP_SD      <- "SD_uniform_root_area"
SHIP_SS      <- "SS_flat"
FORMS        <- expand.grid(sd = SD_FORMS, ss = SS_FORMS, stringsAsFactors = FALSE)
MIN_BLOCK_N  <- 10L   # min structures for a fine block to enter the nbhd model

message(sprintf("[cv] rD=%d cores=%d folds=%d cv_block=%dm nbhd_block=%dm",
                rD, n_cores, N_FOLDS, CV_BLOCK, NBHD_BLOCK))

# ----- Load + terrain join ----------------------------------------------------
sweep <- read_csv(file.path(data_dir, "sweep", "sweep_all.csv"),
                  show_col_types = FALSE) |>
  filter(r_D == rD, !is.na(DAMAGE))

join_fire <- function(fire) {
  kp <- file.path(data_dir, fire, "covariates",
                  sprintf("%s_burned_struc_model_inputs.csv", fire))
  kenny <- read_csv(kp, show_col_types = FALSE)
  ref <- sweep |> filter(fire == !!fire)
  ssdd_sf  <- st_as_sf(ref,   coords = c("cent_x", "cent_y"), crs = EPSG, remove = FALSE)
  kenny_sf <- st_as_sf(kenny, coords = c("utm_x",  "utm_y"),  crs = EPSG)
  nn  <- st_nearest_feature(ssdd_sf, kenny_sf)
  dst <- as.numeric(st_distance(ssdd_sf, kenny_sf[nn, ], by_element = TRUE))
  keep <- dst <= 5
  bind_cols(ref[keep, ], kenny[nn[keep], TERRAIN]) |>
    mutate(destroyed = as.integer(DAMAGE == "Destroyed (>50%)"))
}
dat <- bind_rows(lapply(FIRES, join_fire))
message(sprintf("[cv] %d structures after terrain join (%s)", nrow(dat),
                paste(sprintf("%s=%d", FIRES,
                       vapply(FIRES, \(f) sum(dat$fire == f), integer(1))),
                      collapse = ", ")))

# ----- Coarse spatial fold ids (per fire, shared index) -----------------------
set.seed(7291)
dat$fold <- NA_integer_
for (f in FIRES) {
  idx <- which(dat$fire == f)
  sfd <- st_as_sf(dat[idx, ], coords = c("cent_x", "cent_y"), crs = EPSG)
  fo  <- cv_spatial(x = sfd, k = N_FOLDS, size = CV_BLOCK, selection = "random",
                    iteration = 50, seed = 7291, progress = FALSE,
                    plot = FALSE, report = FALSE)
  dat$fold[idx] <- fo$folds_ids
}

# ----- Fine neighbourhood blocks (per fire), inheriting the coarse fold -------
dat$nblock <- sprintf("%s_%d_%d", dat$fire,
                      as.integer(dat$cent_x %/% NBHD_BLOCK),
                      as.integer(dat$cent_y %/% NBHD_BLOCK))
mode_int <- function(x) as.integer(names(sort(table(x), decreasing = TRUE))[1])
blocks <- dat |>
  group_by(fire, nblock) |>
  summarise(n = n(), rate = mean(destroyed), fold = mode_int(fold),
            across(all_of(c(SD_FORMS, SS_FORMS, TERRAIN)), mean),
            .groups = "drop") |>
  filter(n >= MIN_BLOCK_N)
message(sprintf("[cv] %d fine blocks (>=%d struc) for the neighbourhood model",
                nrow(blocks), MIN_BLOCK_N))

# ----- Scoring helpers --------------------------------------------------------
safe_auc <- function(truth, prob) {
  if (length(unique(truth)) < 2L) return(NA_real_)
  tb <- tibble(t = factor(truth, levels = c(0L, 1L)), p = prob)
  as.numeric(roc_auc_vec(tb$t, tb$p, event_level = "second"))
}
safe_logloss <- function(truth, prob) {
  if (length(unique(truth)) < 2L) return(NA_real_)
  tb <- tibble(t = factor(truth, levels = c(0L, 1L)), p = prob)
  as.numeric(mn_log_loss_vec(tb$t, tb$p, event_level = "second"))
}
concordance <- function(pred, actual) {           # AUC-scale rank agreement
  if (length(pred) < 3L) return(NA_real_)
  dp <- sign(outer(pred, pred, "-")); da <- sign(outer(actual, actual, "-"))
  comp <- da != 0
  if (!any(comp)) return(NA_real_)
  (sum(dp == da & comp) + 0.5 * sum(dp == 0 & comp)) / sum(comp)
}
rmse <- function(pred, actual) sqrt(mean((pred - actual)^2))

fit_structure <- function(feats, train, test) {
  need <- c("destroyed", feats)
  train <- train[complete.cases(train[, need]), ]; test <- test[complete.cases(test[, need]), ]
  if (nrow(train) < 20L || nrow(test) < 1L) return(NULL)
  rf <- ranger(x = train[, feats, drop = FALSE],
               y = factor(train$destroyed, levels = c(0L, 1L)),
               probability = TRUE, num.trees = 500, seed = 7291, num.threads = 1L)
  p <- predict(rf, data = test[, feats, drop = FALSE])$predictions[, "1"]
  data.frame(fire = test$fire, truth = test$destroyed, pred = p)
}
fit_neighbourhood <- function(feats, train, test) {
  need <- c("rate", feats)
  train <- train[complete.cases(train[, need]), ]; test <- test[complete.cases(test[, need]), ]
  if (nrow(train) < 15L || nrow(test) < 3L) return(NULL)
  rf <- ranger(x = train[, feats, drop = FALSE], y = train$rate,
               num.trees = 500, seed = 7291, num.threads = 1L)
  p <- predict(rf, data = test[, feats, drop = FALSE])$predictions
  data.frame(fire = test$fire, actual = test$rate, pred = p)
}

# ----- Build the per-fit work list -------------------------------------------
# Each work item = (scale, cv, form_row, fold). Pooled forms=all 12; within-fire
# secondary uses the shipped form only.
work <- list()
add <- function(...) work[[length(work) + 1L]] <<- list(...)

for (i in seq_len(nrow(FORMS))) {
  sd_f <- FORMS$sd[i]; ss_f <- FORMS$ss[i]
  for (k in seq_len(N_FOLDS)) {
    add(scale = "structure", cv = "pooled", sd = sd_f, ss = ss_f, fold = k)
    add(scale = "neighbourhood", cv = "pooled", sd = sd_f, ss = ss_f, fold = k)
  }
}
for (f in FIRES) for (k in seq_len(N_FOLDS)) {
  add(scale = "structure", cv = "within", sd = SHIP_SD, ss = SHIP_SS, fold = k, only_fire = f)
  add(scale = "neighbourhood", cv = "within", sd = SHIP_SD, ss = SHIP_SS, fold = k, only_fire = f)
}
message(sprintf("[cv] %d fits queued", length(work)))

run_one <- function(w) {
  feats <- c(w$sd, w$ss, TERRAIN)
  restrict <- function(df) if (is.null(w$only_fire)) df else df[df$fire == w$only_fire, ]
  if (w$scale == "structure") {
    d  <- restrict(dat)
    r  <- fit_structure(feats, d[d$fold != w$fold, ], d[d$fold == w$fold, ])
    if (is.null(r) || nrow(r) == 0L) return(NULL)
    grps <- split(r, r$fire); grps[["pooled"]] <- r
    do.call(rbind, lapply(names(grps), function(fr) {
      g <- grps[[fr]]
      data.frame(r_D = rD, scale = "structure", cv = w$cv,
                 sd_form = w$sd, ss_form = w$ss, form = paste(w$sd, w$ss, sep = "|"),
                 fire = fr, fold = w$fold, n = nrow(g),
                 score = safe_auc(g$truth, g$pred), aux = safe_logloss(g$truth, g$pred))
    }))
  } else {
    b  <- restrict(blocks)
    r  <- fit_neighbourhood(feats, b[b$fold != w$fold, ], b[b$fold == w$fold, ])
    if (is.null(r) || nrow(r) == 0L) return(NULL)
    grps <- split(r, r$fire); grps[["pooled"]] <- r
    do.call(rbind, lapply(names(grps), function(fr) {
      g <- grps[[fr]]
      data.frame(r_D = rD, scale = "neighbourhood", cv = w$cv,
                 sd_form = w$sd, ss_form = w$ss, form = paste(w$sd, w$ss, sep = "|"),
                 fire = fr, fold = w$fold, n = nrow(g),
                 score = concordance(g$pred, g$actual), aux = rmse(g$pred, g$actual))
    }))
  }
}

# ----- Run + write ------------------------------------------------------------
t0 <- Sys.time()
res <- mclapply(work, function(w) tryCatch(run_one(w), error = function(e) e),
                mc.cores = n_cores)
errs <- vapply(res, \(x) inherits(x, "error"), logical(1))
if (any(errs)) stop(sprintf("%d fits errored; first: %s", sum(errs),
                            conditionMessage(res[[which(errs)[1]]])))
results <- bind_rows(res[!vapply(res, is.null, logical(1))])
message(sprintf("[cv] %d result rows in %.1fs", nrow(results),
                as.numeric(Sys.time() - t0, units = "secs")))

run_dir <- file.path(out_dir, sprintf("rD%d", rD))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(results, file.path(run_dir, sprintf("rD%d_cv.csv", rD)))
message(sprintf("[cv] wrote %s", file.path(run_dir, sprintf("rD%d_cv.csv", rD))))
