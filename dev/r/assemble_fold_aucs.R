## assemble_fold_aucs.R
##
## Combine the four per-model *_results.csv files into a single
## fold_aucs_rd200_rs50.csv with a "model" column, matching the schema the
## qmd notebook expects.
##
## Run from repo root: Rscript dev/r/assemble_fold_aucs.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

DIR <- "dev/r"

sources <- tribble(
  ~model,          ~path,
  "1_kenny",       file.path(DIR, "benchmark_kenny_rf_results.csv"),
  "2_kenny_ext",   file.path(DIR, "benchmark_kenny_ext_rf_results.csv"),
  "3_ssdd_3met",   file.path(DIR, "benchmark_ssdd_rf_results.csv"),
  "4_ssdd_4met",   file.path(DIR, "baseline_ssdd_rf_results.csv"),
)

fold_aucs <- sources |>
  rowwise() |>
  do({
    read_csv(.$path, show_col_types = FALSE) |>
      mutate(model = .$model, .before = 1) |>
      select(model, setting, fold, auc)
  }) |>
  ungroup()

out <- file.path(DIR, "fold_aucs_rd200_rs50.csv")
write_csv(fold_aucs, out)
message(sprintf("Wrote %s — %d rows across %d models",
                out, nrow(fold_aucs), length(unique(fold_aucs$model))))

fold_aucs |>
  group_by(model, setting) |>
  summarise(mean_auc = mean(auc), .groups = "drop") |>
  arrange(setting, model)
