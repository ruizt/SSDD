## assemble_predictions.R
##
## Combine the four per-model *_predictions.csv files into a single
## predictions_rd200_rs50.csv with a "model" column, matching the schema
## the qmd notebook expects for ROC curves.
##
## Run from repo root: Rscript dev/r/assemble_predictions.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

DIR <- "dev/r"

sources <- tribble(
  ~model,          ~path,
  "1_kenny",       file.path(DIR, "benchmark_kenny_rf_predictions.csv"),
  "2_kenny_ext",   file.path(DIR, "benchmark_kenny_ext_rf_predictions.csv"),
  "3_ssdd_3met",   file.path(DIR, "benchmark_ssdd_rf_predictions.csv"),
  "4_ssdd_4met",   file.path(DIR, "baseline_ssdd_rf_predictions.csv"),
)

preds <- sources |>
  rowwise() |>
  do({
    read_csv(.$path, show_col_types = FALSE) |>
      mutate(model = .$model, .before = 1) |>
      select(model, setting, fold, truth, pred)
  }) |>
  ungroup()

out <- file.path(DIR, "predictions_rd200_rs50.csv")
write_csv(preds, out)
message(sprintf("Wrote %s — %d rows across %d models",
                out, nrow(preds), length(unique(preds$model))))
