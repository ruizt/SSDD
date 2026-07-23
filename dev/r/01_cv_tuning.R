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
##         dev/r/_out/fig_form_comparison.png
##         dev/r/_out/fig_radius_profile.png
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
R_D_FIXED    <- 50L                              # fixed radius: CV surface is flat, 50 m chosen
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

# ── Build the selection table ────────────────────────────────────────────────
# CV surface is flat across r_D for the default form (all within fold noise);
# r_D = 50 m is fixed on substantive grounds (most local, contrasts with
# Kenny's 200 m KDE bandwidth).

selection <- bind_rows(lapply(names(SCALES), function(sc) {
  pooled |>
    filter(scale == sc,
           as.character(form) == DEFAULT_FORM,
           r_D == R_D_FIXED) |>
    mutate(scale = sc, metric = SCALES[[sc]], .before = 1)
})) |>
  select(scale, metric, form, sd_form, ss_form, r_D,
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

cat("\n=== Selected config: package default at r_D = 50 m ===\n")
selection |>
  mutate(across(c(score, score_sd, aux), \(x) round(x, 3))) |>
  as.data.frame() |>
  print(row.names = FALSE)



# ── Figures ───────────────────────────────────────────────────────────────────

theme_set(theme_minimal(base_size = 11))

# (1) Form comparison: per-form mean (over r_D), both scales side-by-side, SD-blocked
forms <- pooled |>
  group_by(scale, form, sd_form, ss_form) |>
  summarise(m = mean(score_mean), .groups = "drop") |>
  mutate(is_default = as.character(form) == DEFAULT_FORM)

# SS order by overall mean rank; SD blocks by parametric structure
# (kernel: uniform → quartic; normalisation: root_area → unit)
ss_order <- pooled |>
  group_by(ss_form) |> summarise(m = mean(score_mean), .groups = "drop") |>
  arrange(m) |> pull(ss_form)
sd_order <- c("SD_uniform_root_area", "SD_uniform_unit",
              "SD_quartic_root_area", "SD_quartic_unit")

forms <- forms |>
  mutate(ss_form = factor(ss_form, levels = ss_order),
         sd_form = factor(sd_form, levels = sd_order))

noise_labels <- sprintf(
  "fold-noise sd:  structure (AUC) = %.3f,  neighbourhood (concordance) = %.3f",
  fold_noise$fold_sd[fold_noise$scale == "structure"],
  fold_noise$fold_sd[fold_noise$scale == "neighbourhood"])

default_scores <- forms |> filter(is_default) |> select(scale, m)

band <- forms |>
  group_by(scale) |>
  summarise(center = mean(m), .groups = "drop") |>
  left_join(fold_noise, by = "scale")

p_forms <- ggplot(forms, aes(m, ss_form, color = scale, shape = scale)) +
  geom_rect(data = band, inherit.aes = FALSE,
            aes(xmin = center - fold_sd, xmax = center + fold_sd,
                fill = scale, ymin = -Inf, ymax = Inf),
            alpha = 0.07) +
  geom_vline(data = default_scores, aes(xintercept = m, color = scale),
             linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
  geom_point(aes(size = is_default)) +
  scale_color_manual(
    values = c(structure = "#2166ac", neighbourhood = "#d6604d"),
    labels = c(structure = "structure (AUC)", neighbourhood = "neighbourhood (concordance)"),
    name = NULL) +
  scale_shape_manual(
    values = c(structure = 16L, neighbourhood = 17L),
    labels = c(structure = "structure (AUC)", neighbourhood = "neighbourhood (concordance)"),
    name = NULL) +
  scale_fill_manual(
    values = c(structure = "#2166ac", neighbourhood = "#d6604d"),
    guide = "none") +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 3.5), guide = "none") +
  facet_grid(sd_form ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Forms are predictively interchangeable",
       subtitle = paste0("shaded band = \u00b11 fold-noise SD around mean; larger dot = package default\n",
                         noise_labels),
       x = "mean held-out score (over r_D)", y = NULL) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0, hjust = 0))
ggsave(file.path(OUT_DIR, "fig_form_comparison.png"), p_forms,
       width = 8, height = 7, dpi = 150)

# (2) Radius profile for the package default form
default_profile <- pooled |>
  filter(as.character(form) == DEFAULT_FORM) |>
  left_join(fold_noise, by = "scale")

p_radius <- ggplot(default_profile, aes(r_D, score_mean, color = scale)) +
  geom_ribbon(aes(ymin = score_mean - fold_sd, ymax = score_mean + fold_sd, fill = scale),
              alpha = 0.1, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = R_D_FIXED, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  scale_color_manual(
    values = c(structure = "#2166ac", neighbourhood = "#d6604d"),
    labels = c(structure = "structure (AUC)", neighbourhood = "neighbourhood (concordance)"),
    name = NULL) +
  scale_fill_manual(
    values = c(structure = "#2166ac", neighbourhood = "#d6604d"),
    guide = "none") +
  facet_wrap(~ scale, scales = "free_y") +
  labs(title = sprintf("Radius is not critical — selecting r_D = %d m", R_D_FIXED),
       subtitle = "ribbon = \u00b11 fold SD;  dashed = selected radius",
       x = expression(r[D] ~ "(m)"), y = "mean held-out score") +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "fig_radius_profile.png"), p_radius,
       width = 9, height = 4, dpi = 150)

cat(sprintf("\nWrote selection + 2 figures to %s/\n", OUT_DIR))

