#!/usr/bin/env Rscript
## collect.R — assemble per-job CV CSVs into one results table, plus a
## per-(r_D, r_S, setting) summary and a per-setting "optimum" table.
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
##   OUT_OPTIMA    best (r_D, r_S) per setting (default <RAW_DIR>/sweep_optima.csv)

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

# ----- Per-setting summary ---------------------------------------------------

summary_tbl <- results |>
  group_by(r_D, r_S, setting) |>
  summarise(
    n_folds       = dplyr::n(),
    auc_mean      = mean(auc,      na.rm = TRUE),
    auc_sd        = sd  (auc,      na.rm = TRUE),
    log_loss_mean = mean(log_loss, na.rm = TRUE),
    brier_mean    = mean(brier,    na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(setting, r_D, r_S)

write_csv(summary_tbl, OUT_SUMMARY)
message(sprintf("Wrote %s", OUT_SUMMARY))

# ----- Per-setting optimum (best held-out AUC) -------------------------------

optima <- summary_tbl |>
  group_by(setting) |>
  slice_max(auc_mean, n = 1L, with_ties = FALSE) |>
  ungroup() |>
  arrange(setting)

write_csv(optima, OUT_OPTIMA)
message(sprintf("Wrote %s", OUT_OPTIMA))

# ----- Print to console ------------------------------------------------------

cat("\n=== Per-setting optimum (best mean held-out AUC across r_D, r_S) ===\n")
print(optima, n = Inf, width = Inf)

cat("\n=== Full summary ===\n")
print(summary_tbl, n = Inf, width = Inf)
