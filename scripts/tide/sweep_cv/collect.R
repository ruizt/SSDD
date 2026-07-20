#!/usr/bin/env Rscript
## collect.R — assemble per-job CV CSVs into one results table, plus a
## per-cell summary and a top-forms "optimum" table.
##
## Run after fetch.sh has pulled per-job outputs to RAW_DIR.
##
## Usage (from the repo root, after fetch):
##   Rscript scripts/tide/sweep_cv/collect.R
##
## Optional env overrides:
##   RAW_DIR       directory holding per-job subdirs (default _data/processed/sweep_cv)
##   OUT_RESULTS   long-format CSV with every fold's score (default <RAW_DIR>/sweep_results.csv)
##   OUT_SUMMARY   mean-per-setting summary (default <RAW_DIR>/sweep_summary.csv)
##   OUT_OPTIMA    top forms per scale (default <RAW_DIR>/sweep_optima.csv)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

RAW_DIR     <- Sys.getenv("RAW_DIR",     "_data/processed/sweep_cv")
OUT_RESULTS <- Sys.getenv("OUT_RESULTS", file.path(RAW_DIR, "sweep_results.csv"))
OUT_SUMMARY <- Sys.getenv("OUT_SUMMARY", file.path(RAW_DIR, "sweep_summary.csv"))
OUT_OPTIMA  <- Sys.getenv("OUT_OPTIMA",  file.path(RAW_DIR, "sweep_optima.csv"))

# ----- Find per-job CSVs ------------------------------------------------------

if (!dir.exists(RAW_DIR)) stop("RAW_DIR not found: ", RAW_DIR)
files <- list.files(RAW_DIR, pattern = "_cv\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0L) stop("No per-job *_cv.csv files under ", RAW_DIR)

message(sprintf("Found %d per-job CV CSVs", length(files)))

# ----- Concatenate ------------------------------------------------------------

results <- map_dfr(files, read_csv, show_col_types = FALSE)
message(sprintf("Concatenated to %d rows", nrow(results)))

dir.create(dirname(OUT_RESULTS), recursive = TRUE, showWarnings = FALSE)
write_csv(results, OUT_RESULTS)
message(sprintf("Wrote %s", OUT_RESULTS))

# ----- Per-cell summary (mean over folds) ------------------------------------
# One row per (r_D, scale, cv, form, fire). score = AUC (structure) /
# concordance (neighbourhood); aux = log_loss / rmse.

summary_tbl <- results |>
  group_by(r_D, scale, cv, sd_form, ss_form, form, fire) |>
  summarise(
    n_folds    = dplyr::n(),
    score_mean = mean(score, na.rm = TRUE),
    score_sd   = sd  (score, na.rm = TRUE),
    aux_mean   = mean(aux,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(scale, cv, fire, form, r_D)

write_csv(summary_tbl, OUT_SUMMARY)
message(sprintf("Wrote %s", OUT_SUMMARY))

# ----- Optima: best (r_D, form) per scale, on the pooled held-out set --------

optima <- summary_tbl |>
  filter(cv == "pooled", fire == "pooled") |>
  group_by(scale) |>
  slice_max(score_mean, n = 5L, with_ties = FALSE) |>
  ungroup() |>
  arrange(scale, desc(score_mean))

write_csv(optima, OUT_OPTIMA)
message(sprintf("Wrote %s", OUT_OPTIMA))

# ----- Print to console ------------------------------------------------------

cat("\n=== Top (r_D, form) per scale on the pooled held-out set ===\n")
print(optima, n = Inf, width = Inf)

cat("\n=== Shipped form (SD_uniform_root_area | SS_flat), by scale/fire/r_D ===\n")
summary_tbl |>
  filter(cv == "pooled", form == "SD_uniform_root_area|SS_flat") |>
  arrange(scale, fire, r_D) |>
  print(n = Inf, width = Inf)
