## 01_cv_tuning.R
##
## Post-process the Tide sweep_cv results into a hyperparameter-tuning summary.
##
## The CV sweep (scripts/tide/sweep_cv) fit, for every r_D and every metric
## form-setting (4 SD × 3 SS = 12), a {SD, SS} + terrain random forest under one
## pooled spatial-block CV, at two scales (structure classifier, neighbourhood
## aggregate-then-fit). This script treats that as a tuning grid and answers:
##
##   (a) how does held-out skill move across form and radius?   -> tables + figures
##   (b) what is the tuned skill of
##         - the PACKAGE DEFAULT form  at its best radius
##         - the OPTIMAL form          at its best radius
##       at each scale, with a per-fire breakdown.
##
## "Best radius" = the r_D maximising pooled held-out score for that form.
## "Optimal form" = the (form, r_D) cell maximising pooled held-out score.
## Score = AUC (structure) or concordance (neighbourhood); both higher-is-better.
##
## The two selected configs are written to dev/r/_out/tuning_selection.csv, which
## 02_benchmark.R consumes to ship the tuned models into the Kenny head-to-head.
##
## Input : _data/processed/sweep_cv/sweep_summary.csv  (mean over folds)
##         _data/processed/sweep_cv/sweep_results.csv   (per-fold, for noise band)
## Output: dev/r/_out/tuning_selection.csv
##         dev/r/_out/fig_form_radius_heatmap.png
##         dev/r/_out/fig_radius_profile.png
##         dev/r/_out/fig_form_flatness.png
##         dev/r/_out/fig_fire_disaggregation.png
##
## Run from repo root:  Rscript dev/r/01_cv_tuning.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ── Config ────────────────────────────────────────────────────────────────────

CV_DIR       <- "_data/processed/sweep_cv"
OUT_DIR      <- "dev/r/_out"
DEFAULT_FORM <- "SD_uniform_root_area|SS_flat"   # the shipped package default
SCALES       <- c(structure = "AUC", neighbourhood = "concordance")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Load ──────────────────────────────────────────────────────────────────────

summary_tbl <- read_csv(file.path(CV_DIR, "sweep_summary.csv"), show_col_types = FALSE)
results_tbl <- read_csv(file.path(CV_DIR, "sweep_results.csv"),  show_col_types = FALSE)

stopifnot(DEFAULT_FORM %in% summary_tbl$form)

# Pooled held-out surface: one score per (scale, form, r_D).
pooled <- summary_tbl |>
  filter(cv == "pooled", fire == "pooled") |>
  select(scale, form, sd_form, ss_form, r_D, score_mean, score_sd, aux_mean)

# Typical fold-to-fold noise, per scale — the yardstick for "are forms different".
fold_noise <- results_tbl |>
  filter(cv == "pooled", fire == "pooled") |>
  group_by(scale) |>
  summarise(fold_sd = sd(score, na.rm = TRUE), .groups = "drop")

# ── Selection helpers ─────────────────────────────────────────────────────────

# Best radius for a fixed form, at a given scale.
best_radius_for_form <- function(scale_name, form_name) {
  pooled |>
    filter(scale == scale_name, form == form_name) |>
    slice_max(score_mean, n = 1, with_ties = FALSE)
}

# Best (form, r_D) cell at a given scale.
best_cell <- function(scale_name) {
  pooled |>
    filter(scale == scale_name) |>
    slice_max(score_mean, n = 1, with_ties = FALSE)
}

# Per-fire scores at a chosen (scale, form, r_D).
per_fire_at <- function(scale_name, form_name, r_d) {
  summary_tbl |>
    filter(cv == "pooled", scale == scale_name, form == form_name, r_D == r_d,
           fire != "pooled") |>
    select(fire, score_mean) |>
    arrange(fire)
}

# ── Build the selection table ────────────────────────────────────────────────

selection <- bind_rows(lapply(names(SCALES), function(sc) {
  d <- best_radius_for_form(sc, DEFAULT_FORM) |> mutate(config = "package_default")
  o <- best_cell(sc)                          |> mutate(config = "optimal_form")
  bind_rows(d, o) |>
    mutate(scale = sc, metric = SCALES[[sc]], .before = 1)
}))

selection <- selection |>
  select(scale, metric, config, form, sd_form, ss_form, r_D,
         score = score_mean, score_sd, aux = aux_mean)

write_csv(selection, file.path(OUT_DIR, "tuning_selection.csv"))

# ── Console report ────────────────────────────────────────────────────────────

cat("\n=================================================================\n")
cat("  SSDD hyperparameter tuning — pooled held-out CV surface\n")
cat("=================================================================\n")

for (sc in names(SCALES)) {
  metric <- SCALES[[sc]]
  noise  <- fold_noise$fold_sd[fold_noise$scale == sc]
  surf   <- pooled |> filter(scale == sc)
  cat(sprintf("\n--- %s scale (%s) ---\n", toupper(sc), metric))
  cat(sprintf("  form spread:  best %.3f  worst %.3f  range %.3f   (fold noise sd %.3f)\n",
              max(surf$score_mean), min(surf$score_mean),
              max(surf$score_mean) - min(surf$score_mean), noise))
  cat(sprintf("  radius spread (default form): %.3f - %.3f across r_D %d-%d\n",
              min(surf$score_mean[surf$form == DEFAULT_FORM]),
              max(surf$score_mean[surf$form == DEFAULT_FORM]),
              min(surf$r_D), max(surf$r_D)))
}

cat("\n=== Tuned selection (each config at its best radius) ===\n")
selection |>
  mutate(across(c(score, score_sd, aux), \(x) round(x, 3))) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\n=== Per-fire held-out score at each selected config ===\n")
for (i in seq_len(nrow(selection))) {
  r <- selection[i, ]
  pf <- per_fire_at(r$scale, r$form, r$r_D)
  cat(sprintf("  [%s | %s] %s @ r_D=%d  (pooled %.3f):  %s\n",
              r$scale, r$config, r$form, r$r_D, r$score,
              paste(sprintf("%s=%.3f", pf$fire, round(pf$score_mean, 3)), collapse = "  ")))
}

# ── Figures ───────────────────────────────────────────────────────────────────

theme_set(theme_minimal(base_size = 11))

# order forms by overall mean score (average across scales) for stable y-axis
form_order <- pooled |>
  group_by(form) |> summarise(m = mean(score_mean), .groups = "drop") |>
  arrange(m) |> pull(form)

pooled <- pooled |> mutate(form = factor(form, levels = form_order))

# (1) Heatmap: form × r_D, faceted by scale
p_heat <- ggplot(pooled, aes(factor(r_D), form, fill = score_mean)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", score_mean)), size = 2.6) +
  scale_fill_viridis_c(option = "mako", name = "held-out\nscore") +
  facet_wrap(~ scale, scales = "free_x") +
  labs(title = "Held-out skill across metric form and radius",
       subtitle = "pooled spatial-block CV; structure = AUC, neighbourhood = concordance",
       x = expression(r[D] ~ "(m)"), y = NULL) +
  theme(panel.grid = element_blank())
ggsave(file.path(OUT_DIR, "fig_form_radius_heatmap.png"), p_heat,
       width = 11, height = 5.5, dpi = 150)

# (2) Radius profile: default vs optimal-form line per scale
profile <- pooled |>
  mutate(highlight = case_when(
    as.character(form) == DEFAULT_FORM ~ "package default",
    TRUE ~ "other forms"))
opt_forms <- selection |> filter(config == "optimal_form") |>
  distinct(scale, form) |> mutate(opt = TRUE)
profile <- profile |>
  left_join(opt_forms, by = c("scale", "form")) |>
  mutate(highlight = if_else(!is.na(opt), "optimal form", highlight))

p_prof <- ggplot(profile, aes(r_D, score_mean, group = form)) +
  geom_line(data = \(d) filter(d, highlight == "other forms"),
            color = "grey80", linewidth = 0.4) +
  geom_line(data = \(d) filter(d, highlight != "other forms"),
            aes(color = highlight), linewidth = 1) +
  geom_point(data = \(d) filter(d, highlight != "other forms"),
             aes(color = highlight), size = 2) +
  scale_color_manual(values = c("package default" = "#1b7837",
                                "optimal form" = "#762a83"), name = NULL) +
  facet_wrap(~ scale, scales = "free_y") +
  labs(title = "Radius profile: default and optimal forms vs all others",
       subtitle = "grey = the other 10 form-settings (near-identical)",
       x = expression(r[D] ~ "(m)"), y = "pooled held-out score") +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "fig_radius_profile.png"), p_prof,
       width = 10, height = 5, dpi = 150)

# (3) Form flatness: per-form mean (over r_D) with fold-noise band, default marked
flat <- pooled |>
  group_by(scale, form) |>
  summarise(m = mean(score_mean), .groups = "drop") |>
  left_join(fold_noise, by = "scale") |>
  mutate(is_default = as.character(form) == DEFAULT_FORM)
band <- flat |> group_by(scale) |>
  summarise(center = mean(m), noise = first(fold_sd), .groups = "drop")

p_flat <- ggplot(flat, aes(m, form)) +
  geom_rect(data = band, inherit.aes = FALSE,
            aes(xmin = center - noise, xmax = center + noise,
                ymin = -Inf, ymax = Inf), alpha = 0.12, fill = "grey50") +
  geom_point(aes(color = is_default, size = is_default)) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#1b7837"),
                     labels = c("form", "package default"), name = NULL) +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 3.5), guide = "none") +
  facet_wrap(~ scale, scales = "free_x") +
  labs(title = "Forms are predictively interchangeable",
       subtitle = "grey band = ±1 fold-to-fold sd around the mean form; all forms fall inside it",
       x = "mean held-out score (over r_D)", y = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "fig_form_flatness.png"), p_flat,
       width = 10, height = 5, dpi = 150)

# (4) Fire disaggregation at the selected configs
fire_rows <- bind_rows(lapply(seq_len(nrow(selection)), function(i) {
  r <- selection[i, ]
  per_fire_at(r$scale, r$form, r$r_D) |>
    mutate(scale = r$scale, config = r$config, r_D = r$r_D)
}))

p_fire <- ggplot(fire_rows, aes(reorder(fire, score_mean), score_mean, fill = config)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = sprintf("%.2f", score_mean)),
            position = position_dodge(0.8), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("package_default" = "#1b7837",
                               "optimal_form" = "#762a83"), name = NULL) +
  facet_wrap(~ scale) +
  labs(title = "Per-fire held-out skill at the tuned configs",
       subtitle = "dashed line = chance (0.5); the pooled model is carried by Palisades",
       x = NULL, y = "held-out score") +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "fig_fire_disaggregation.png"), p_fire,
       width = 9, height = 5, dpi = 150)

cat(sprintf("\nWrote selection + 4 figures to %s/\n", OUT_DIR))
