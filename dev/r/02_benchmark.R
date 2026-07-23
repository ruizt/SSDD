## 02_benchmark.R
##
## Downstream head-to-head: does the tuned SSDD package beat the existing
## (no-SSDD) building covariates on held-out damage prediction?
##
## Ships the tuned configs chosen by 01_cv_tuning.R (dev/r/_out/tuning_selection.csv,
## structure scale) and fits three structure-level RF classifiers on identical
## spatial-block folds and identical observations:
##
##   kenny_ext     Kenny's building covariates + terrain  (the no-SSDD reference)
##                   build_dens, area_ext_build_2, distance_to_nearest_building,
##                   elevation, slope, aspect
##   ssdd_default  tuned PACKAGE DEFAULT form + terrain
##                   e.g. SD_uniform_root_area, SS_flat  (at its best r_D) + terrain
##   ssdd_optimal  tuned OPTIMAL form + terrain
##                   e.g. SD_uniform_unit, SS_gauss      (at its best r_D) + terrain
##
## Only Eaton and Palisades carry Kenny building covariates (Mountain's covariate
## file is terrain-only, from the DEM), so the comparison is those two fires.
## The CV machinery mirrors scripts/tide/sweep_cv/compute.R exactly (same folds,
## same label, same ranger settings), so ssdd_* here reproduce the sweep_cv AUCs
## as an internal check.
##
## NOTE on features: kenny_ext's distance_to_nearest_building is LITERAL metres
## (higher = more isolated); the package's SS_flat is its monotone inverse,
## 1/(d_nn+eps). For a rank-split RF these carry the same separation information,
## so the SSDD novelty over kenny_ext is really the SD density kernel (+ the
## orientation option in the optimal form), not the nearest-neighbour term.
##
## Input : _data/processed/sweep/sweep_all.csv
##         _data/processed/<fire>/covariates/<fire>_burned_struc_model_inputs.csv
##         dev/r/_out/tuning_selection.csv        (from 01_cv_tuning.R)
## Output: dev/r/_out/benchmark_results.csv       (per model × fold)
##         dev/r/_out/benchmark_predictions.csv   (per model × held-out row)
##         dev/r/_out/fig_benchmark_roc.png
##         dev/r/_out/fig_benchmark_auc.png
##
## The 30 model×fold fits run in parallel (parallel::mclapply, all cores by
## default; set SSDD_CORES to cap). ranger stays single-threaded in each fork.
##
## Run from repo root:  Rscript dev/r/02_benchmark.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(sf)
  library(blockCV)
  library(ranger)
  library(yardstick)
  library(ggplot2)
  library(patchwork)
  library(parallel)
})

# ── Config ────────────────────────────────────────────────────────────────────

PROCESSED <- "_data/processed"
SWEEP_ALL <- file.path(PROCESSED, "sweep", "sweep_all.csv")
SELECTION <- "dev/r/_out/tuning_selection.csv"
OUT_DIR   <- "dev/r/_out"
FIRES     <- c("eaton", "palisades")        # Kenny building covariates exist here only
EPSG      <- 32611L
N_FOLDS   <- 10L
CV_BLOCK  <- 1500                            # matches sweep_cv coarse fold block
SEED      <- 7291L
TERRAIN   <- c("elevation", "slope", "aspect")
KENNY_BLD <- c("build_dens", "area_ext_build_2", "distance_to_nearest_building")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

# ── 1. Tuned configs from 01_cv_tuning.R (structure scale) ────────────────────

sel <- read_csv(SELECTION, show_col_types = FALSE) |> filter(scale == "structure")
def_row <- sel |> filter(config == "package_default")
opt_row <- sel |> filter(config == "optimal_form")
stopifnot(nrow(def_row) == 1L, nrow(opt_row) == 1L)

message(sprintf("[bench] default: %s @ r_D=%d   |   optimal: %s @ r_D=%d",
                def_row$form, def_row$r_D, opt_row$form, opt_row$r_D))

# ── 2. SSDD feature slices at each tuned r_D ──────────────────────────────────

sweep <- read_csv(SWEEP_ALL, show_col_types = FALSE) |>
  filter(fire %in% FIRES, !is.na(DAMAGE))

def_feats <- sweep |> filter(r_D == def_row$r_D) |>
  transmute(fire, ssdd_id,
            SD_def = .data[[def_row$sd_form]], SS_def = .data[[def_row$ss_form]])
opt_feats <- sweep |> filter(r_D == opt_row$r_D) |>
  transmute(fire, ssdd_id,
            SD_opt = .data[[opt_row$sd_form]], SS_opt = .data[[opt_row$ss_form]])

# Unique buildings (label + coords) from a single r_D slice.
base <- sweep |> filter(r_D == def_row$r_D) |>
  distinct(fire, ssdd_id, cent_x, cent_y, DAMAGE)

# ── 3. NN-join Kenny building covariates + terrain (<=5 m), per fire ──────────

join_fire <- function(fire) {
  kp <- file.path(PROCESSED, fire, "covariates",
                  sprintf("%s_burned_struc_model_inputs.csv", fire))
  kenny <- read_csv(kp, show_col_types = FALSE)
  need  <- c(KENNY_BLD, TERRAIN)
  ref <- base |> filter(fire == !!fire)
  ssdd_sf  <- st_as_sf(ref,   coords = c("cent_x", "cent_y"), crs = EPSG, remove = FALSE)
  kenny_sf <- st_as_sf(kenny, coords = c("utm_x",  "utm_y"),  crs = EPSG)
  nn  <- st_nearest_feature(ssdd_sf, kenny_sf)
  dst <- as.numeric(st_distance(ssdd_sf, kenny_sf[nn, ], by_element = TRUE))
  keep <- dst <= 5
  bind_cols(ref[keep, ], kenny[nn[keep], need]) |>
    mutate(destroyed = as.integer(DAMAGE == "Destroyed (>50%)"))
}

dat <- bind_rows(lapply(FIRES, join_fire)) |>
  left_join(def_feats, by = c("fire", "ssdd_id")) |>
  left_join(opt_feats, by = c("fire", "ssdd_id"))

message(sprintf("[bench] %d structures after Kenny join (%s); %d destroyed",
                nrow(dat),
                paste(sprintf("%s=%d", FIRES,
                       vapply(FIRES, \(f) sum(dat$fire == f), integer(1))),
                      collapse = ", "),
                sum(dat$destroyed)))

# ── 4. Shared spatial-block folds (per fire, common index) ────────────────────

dat$fold <- NA_integer_
for (f in FIRES) {
  idx <- which(dat$fire == f)
  sfd <- st_as_sf(dat[idx, ], coords = c("cent_x", "cent_y"), crs = EPSG)
  fo  <- cv_spatial(x = sfd, k = N_FOLDS, size = CV_BLOCK, selection = "random",
                    iteration = 50, seed = SEED, progress = FALSE,
                    plot = FALSE, report = FALSE)
  dat$fold[idx] <- fo$folds_ids
}

# ── 5. Model definitions + CV fit ─────────────────────────────────────────────

MODELS <- list(
  kenny_ext    = c(KENNY_BLD, TERRAIN),
  ssdd_default = c("SD_def", "SS_def", TERRAIN),
  ssdd_optimal = c("SD_opt", "SS_opt", TERRAIN)
)

fit_predict <- function(feats, train, test) {
  need  <- c("destroyed", feats)
  train <- train[complete.cases(train[, need]), ]
  test  <- test [complete.cases(test [, need]), ]
  if (nrow(train) < 20L || nrow(test) < 1L) return(NULL)
  rf <- ranger(x = train[, feats, drop = FALSE],
               y = factor(train$destroyed, levels = c(0L, 1L)),
               probability = TRUE, num.trees = 500, seed = SEED, num.threads = 1L)
  p <- predict(rf, data = test[, feats, drop = FALSE])$predictions[, "1"]
  data.frame(fire = test$fire, truth = test$destroyed, pred = p)
}

# The model × fold fits are independent; each ranger is seeded, so the result
# is order-independent and reproducible. Fan them out across cores (ranger stays
# single-threaded inside each fork). Override with SSDD_CORES=1 to force serial.
N_CORES <- as.integer(Sys.getenv("SSDD_CORES", parallel::detectCores()))
tasks   <- expand.grid(model = names(MODELS), fold = seq_len(N_FOLDS),
                       stringsAsFactors = FALSE)
message(sprintf("[bench] %d fits (%d models × %d folds) on %d cores",
                nrow(tasks), length(MODELS), N_FOLDS, N_CORES))

preds <- bind_rows(mclapply(seq_len(nrow(tasks)), function(i) {
  m <- tasks$model[i]; k <- tasks$fold[i]
  out <- fit_predict(MODELS[[m]], dat[dat$fold != k, ], dat[dat$fold == k, ])
  if (is.null(out)) return(NULL)
  out$model <- m; out$fold <- k; out
}, mc.cores = N_CORES))

write_csv(preds, file.path(OUT_DIR, "benchmark_predictions.csv"))

# ── 6. Scores ─────────────────────────────────────────────────────────────────

auc_of  <- function(truth, pred) {
  if (length(unique(truth)) < 2L) return(NA_real_)
  as.numeric(roc_auc_vec(factor(truth, levels = c(0L, 1L)), pred, event_level = "second"))
}
ll_of <- function(truth, pred) {
  if (length(unique(truth)) < 2L) return(NA_real_)
  as.numeric(mn_log_loss_vec(factor(truth, levels = c(0L, 1L)), pred, event_level = "second"))
}

# per (model, fold)
per_fold <- preds |>
  group_by(model, fold) |>
  summarise(auc = auc_of(truth, pred), log_loss = ll_of(truth, pred),
            n = n(), .groups = "drop")

# per (model, fire) pooled over folds
per_fire <- preds |>
  group_by(model, fire) |>
  summarise(auc = auc_of(truth, pred), .groups = "drop")

# per model: mean fold AUC + pooled AUC
per_model <- per_fold |>
  group_by(model) |>
  summarise(auc_mean = mean(auc, na.rm = TRUE),
            auc_sd   = sd(auc, na.rm = TRUE),
            ll_mean  = mean(log_loss, na.rm = TRUE), .groups = "drop")
pooled_auc <- preds |> group_by(model) |>
  summarise(auc_pooled = auc_of(truth, pred), .groups = "drop")
per_model <- left_join(per_model, pooled_auc, by = "model")

results <- per_fold |> mutate(scope = "fold") |>
  bind_rows(per_fire |> transmute(model, fold = NA_integer_, auc,
                                  scope = paste0("fire:", fire)))
write_csv(results, file.path(OUT_DIR, "benchmark_results.csv"))

# Paired (same-fold) deltas vs the Kenny reference.
wide <- per_fold |> select(model, fold, auc) |>
  pivot_wider(names_from = model, values_from = auc)
paired <- tibble(
  contrast = c("ssdd_default - kenny_ext", "ssdd_optimal - kenny_ext"),
  mean_delta = c(mean(wide$ssdd_default - wide$kenny_ext, na.rm = TRUE),
                 mean(wide$ssdd_optimal - wide$kenny_ext, na.rm = TRUE)),
  sd_delta   = c(sd(wide$ssdd_default - wide$kenny_ext, na.rm = TRUE),
                 sd(wide$ssdd_optimal - wide$kenny_ext, na.rm = TRUE)))

# ── 7. Console report ─────────────────────────────────────────────────────────

cat("\n=================================================================\n")
cat("  Downstream benchmark — tuned SSDD vs extended-Kenny (no SSDD)\n")
cat("  Fires: Eaton + Palisades | ", N_FOLDS, "spatial-block folds\n")
cat("=================================================================\n\n")
cat("=== Per model (held-out AUC) ===\n")
per_model |> mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  as.data.frame() |> print(row.names = FALSE)

cat("\n=== Per model × fire (pooled over folds) ===\n")
per_fire |> pivot_wider(names_from = fire, values_from = auc) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  as.data.frame() |> print(row.names = FALSE)

cat("\n=== Paired same-fold deltas vs kenny_ext ===\n")
paired |> mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  as.data.frame() |> print(row.names = FALSE)

# ── 8. Figures ────────────────────────────────────────────────────────────────

theme_set(theme_minimal(base_size = 11))
pal <- c(kenny_ext = "#b2182b", ssdd_default = "#1b7837", ssdd_optimal = "#762a83")

# (a) pooled ROC curves
roc_df <- preds |>
  mutate(truth_f = factor(truth, levels = c(0L, 1L))) |>
  group_by(model) |>
  roc_curve(truth = truth_f, pred, event_level = "second")
lab <- per_model |>
  mutate(label = sprintf("%s (AUC %.3f)", model, auc_pooled)) |>
  select(model, label)
roc_df <- left_join(roc_df, lab, by = "model")

p_roc <- ggplot(roc_df, aes(1 - specificity, sensitivity, color = label)) +
  geom_abline(linetype = "dashed", color = "grey60") +
  geom_path(linewidth = 1) +
  scale_color_manual(values = setNames(pal[lab$model], lab$label), name = NULL) +
  coord_equal() +
  labs(title = "Held-out ROC: tuned SSDD vs extended-Kenny",
       subtitle = "pooled over folds, Eaton + Palisades",
       x = "1 - specificity", y = "sensitivity") +
  theme(legend.position = c(0.62, 0.18))

# (b) per-fold AUC dotplot + per-fire
p_auc <- ggplot(per_fold, aes(reorder(model, auc, mean), auc, color = model)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_jitter(width = 0.12, height = 0, alpha = 0.5, size = 1.6) +
  stat_summary(fun = mean, geom = "point", size = 3.5, shape = 18, color = "black") +
  scale_color_manual(values = pal, guide = "none") +
  coord_flip() +
  labs(title = "Per-fold held-out AUC", subtitle = "black diamond = mean",
       x = NULL, y = "AUC")

ggsave(file.path(OUT_DIR, "fig_benchmark_roc.png"), p_roc, width = 6.5, height = 6, dpi = 150)
ggsave(file.path(OUT_DIR, "fig_benchmark_auc.png"), p_auc, width = 7, height = 4, dpi = 150)

cat(sprintf("\nWrote results, predictions, and 2 figures to %s/\n", OUT_DIR))
