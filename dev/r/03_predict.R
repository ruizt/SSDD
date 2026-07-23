## 03_predict.R
##
## Out-of-sample structure-level damage prediction: fit ssdd_default RF on
## Palisades + Mountain, predict on Eaton.  Produces a point-layer GeoPackage
## of P(destroyed) per Eaton structure — a genuine out-of-fire prediction
## mimicking the deployment workflow (train on seen fires, predict on a new
## community).
##
## Mountain SSDD features are in per-fire sweep directories rather than
## sweep_all.csv, so they are loaded separately.  Terrain for all fires is
## NN-joined (≤5 m) from each fire's burned-structure model-inputs CSV.
##
## Input : _data/processed/sweep/sweep_all.csv            (Eaton + Palisades)
##         _data/processed/sweep/mountain_rD50/           (Mountain at r_D=50)
##         _data/processed/<fire>/covariates/<fire>_burned_struc_model_inputs.csv
##         dev/r/_out/tuning_selection.csv
## Output: dev/r/_out/eaton_predicted_damage.gpkg
##
## Run from repo root:  Rscript dev/r/03_predict.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(ranger)
})

# ── Config ────────────────────────────────────────────────────────────────────

PROCESSED   <- "_data/processed"
SWEEP_ALL   <- file.path(PROCESSED, "sweep", "sweep_all.csv")
SELECTION   <- "dev/r/_out/tuning_selection.csv"
OUT_DIR     <- "dev/r/_out"
TRAIN_FIRES <- c("palisades", "mountain")
PRED_FIRE   <- "eaton"
EPSG        <- 32611L
TERRAIN     <- c("elevation", "slope", "aspect")
SEED        <- 7291L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Tuned config ───────────────────────────────────────────────────────────

def_row <- read_csv(SELECTION, show_col_types = FALSE) |>
  filter(scale == "structure")
stopifnot(nrow(def_row) == 1L)

SD_COL <- def_row$sd_form   # "SD_uniform_root_area"
SS_COL <- def_row$ss_form   # "SS_flat"
R_D    <- def_row$r_D       # 50

message(sprintf("[predict] %s | %s @ r_D=%d", SD_COL, SS_COL, R_D))

# ── 2. SSDD features ──────────────────────────────────────────────────────────

# Eaton + Palisades from sweep_all.csv
sweep_ep <- read_csv(SWEEP_ALL, show_col_types = FALSE) |>
  filter(r_D == R_D, !is.na(DAMAGE)) |>
  transmute(fire, ssdd_id, cent_x, cent_y, DAMAGE,
            SD = .data[[SD_COL]], SS = .data[[SS_COL]])

# Mountain from its per-fire sweep directory
sweep_mtn <- read_csv(
  file.path(PROCESSED, "sweep", sprintf("mountain_rD%d", R_D),
            sprintf("mountain_rD%d_metrics.csv", R_D)),
  show_col_types = FALSE) |>
  filter(!is.na(DAMAGE)) |>
  transmute(fire = "mountain", ssdd_id, cent_x, cent_y, DAMAGE,
            SD = .data[[SD_COL]], SS = .data[[SS_COL]])

feats <- bind_rows(sweep_ep, sweep_mtn)

message(sprintf("[predict] %d structures with DAMAGE labels (%s)",
                nrow(feats),
                paste(sprintf("%s=%d", unique(feats$fire),
                              tabulate(factor(feats$fire, unique(feats$fire)))),
                      collapse = ", ")))

# ── 3. Join terrain per fire (NN ≤5 m) ───────────────────────────────────────

join_terrain <- function(fire_name) {
  kp  <- file.path(PROCESSED, fire_name, "covariates",
                   sprintf("%s_burned_struc_model_inputs.csv", fire_name))
  cov <- read_csv(kp, show_col_types = FALSE)
  ref <- feats |> filter(fire == fire_name)

  ref_sf <- st_as_sf(ref, coords = c("cent_x", "cent_y"), crs = EPSG, remove = FALSE)
  cov_sf <- st_as_sf(cov, coords = c("utm_x",  "utm_y"),  crs = EPSG)

  nn  <- st_nearest_feature(ref_sf, cov_sf)
  dst <- as.numeric(st_distance(ref_sf, cov_sf[nn, ], by_element = TRUE))
  keep <- dst <= 5

  bind_cols(ref[keep, ], cov[nn[keep], TERRAIN]) |>
    mutate(destroyed = as.integer(DAMAGE == "Destroyed (>50%)"))
}

all_fires <- unique(feats$fire)
dat <- bind_rows(lapply(all_fires, join_terrain))

message(sprintf("[predict] after terrain join: %d structures (%s)",
                nrow(dat),
                paste(sprintf("%s=%d", all_fires,
                              vapply(all_fires, \(f) sum(dat$fire == f), integer(1))),
                      collapse = ", ")))

# ── 4. Fit on Palisades + Mountain ───────────────────────────────────────────

FEATS <- c("SD", "SS", TERRAIN)

train <- dat |>
  filter(fire %in% TRAIN_FIRES) |>
  filter(if_all(all_of(c("destroyed", FEATS)), \(x) !is.na(x)))

message(sprintf("[predict] training on %d structures (%d destroyed, %d intact)",
                nrow(train), sum(train$destroyed), sum(!train$destroyed)))

rf <- ranger(x     = train[, FEATS, drop = FALSE],
             y     = factor(train$destroyed, levels = c(0L, 1L)),
             probability = TRUE, num.trees = 500, seed = SEED)

# ── 5. Predict on Eaton ───────────────────────────────────────────────────────

pred_dat <- dat |>
  filter(fire == PRED_FIRE) |>
  filter(if_all(all_of(FEATS), \(x) !is.na(x)))

pred_dat$p_destroyed <- predict(rf, data = pred_dat[, FEATS, drop = FALSE])$predictions[, "1"]

message(sprintf("[predict] predicted %d Eaton structures  (mean P=%.3f, median P=%.3f)",
                nrow(pred_dat),
                mean(pred_dat$p_destroyed),
                median(pred_dat$p_destroyed)))

# ── 6. Write GeoPackage ───────────────────────────────────────────────────────

out_sf <- pred_dat |>
  select(ssdd_id, DAMAGE, destroyed, p_destroyed, cent_x, cent_y) |>
  st_as_sf(coords = c("cent_x", "cent_y"), crs = EPSG)

out_path <- file.path(OUT_DIR, "eaton_predicted_damage.gpkg")
st_write(out_sf, out_path, delete_dsn = TRUE, quiet = TRUE)
message(sprintf("[predict] wrote %s", out_path))
