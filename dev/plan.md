# SSDD Development Plan

Concrete steps for developing the codes that consume the ssdd package
outputs.

## 1. Scrape terrain (slope, aspect, elevation) from USGS 3DEP

R-side workflow:

- Read the per-fire raw-metrics CSV; convert centroids to `sf` points.
- Fetch DEM tiles via `elevatr::get_elev_raster` (z=13, ≈10 m).
- Compute slope and aspect with `terra::terrain`.
- Extract values at each centroid; join to the raw-metrics table by
  `ssdd_id`.

Default source is public 3DEP — uniform across all fires, no
collaborator dependency. Higher-resolution LiDAR DTMs can be substituted
later as a sensitivity check if they become available.

## 2. Build the spatial cross-validation harness in R

Foundational infrastructure that downstream modeling steps consume.

- **Within-fire CV**: `blockCV::spatialBlock` with ~500 m blocks,
  tuned from a variogram on residuals once a baseline model exists.
  Five folds.
- **Cross-fire CV**: leave-one-fire-out (LOFO) splits over the pooled
  dataset.
- **Blend re-fitting**: any blend learned on TRAIN must be re-fit
  within each fold (not on the full dataset). Folds are responsible
  for guarding against leakage.
- **Reported metrics**: AUC, log-loss, Brier per fold; mean ± SD across
  folds; per-fire and pooled.

The harness is written as a small set of reusable R functions in `dev/`
so it can be called from each subsequent step without reimplementation.

## 3. Two-step PLS + RF baseline (per-fire and pooled)

For each spatial-CV fold:

1. Residualize `DAMAGE` against controls on TRAIN.
2. Residualize each raw SSDD metric (`KD_raw`, `BA_raw`, `DP_raw`,
   `OP_raw`) against controls on TRAIN.
3. Fit PLS of damage-residual on metric-residuals → one latent component;
   the loading vector is the learned blend.
4. Apply the blend to TEST → single SSDD score.
5. Fit `ranger` on `(blend_score, controls)` over TRAIN; predict TEST.
6. Score.

Both per-fire and pooled fits are run. Comparing blend weights and
held-out skill between configurations characterizes cross-fire
heterogeneity directly.

## 4. Grouped oblique random forest prototype

Joint estimator with two feature groups:

- **Blend group** — the four raw SSDD metrics. Each split is a
  PLS-derived oblique direction (response-aware).
- **Control group** — all other covariates. Each split is axis-aligned,
  as in standard RF.

Bagging and random feature subsetting apply as usual. Evaluated under the
same spatial-CV harness from step 2, per fire and pooled, with the
step-3 two-step pipeline as the comparison baseline.

Implementation: first check whether `ODRF` (or `RotationForest`) can be
extended with grouped-split logic. Fall back to writing from scratch on
top of `rpart`'s tree-building scaffolding if extension is harder than
clean implementation.

## 5. Migrate R modeling code into an R package

Once the workflow stabilizes, promote the R code from `dev/` to a proper
R package (`Rssdd/`) with `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, and
tests. Lives alongside the Python package; both install independently
from the repo root.

## R modeling stack

| Package | Use |
|---|---|
| `ranger` | Fast axis-aligned random forest |
| `pls` | Partial least squares regression |
| `blockCV` | Spatial block cross-validation folds |
| `sf` | Spatial I/O and geometry |
| `terra` | Raster I/O for terrain |
| `elevatr` | USGS 3DEP DEM download |
| `tidyverse` / `data.table` | Data manipulation |
| `quarto` / `rmarkdown` | Reproducible workflow docs |

Situational: `Renvlp` (envelope estimators), `ODRF` (oblique RF base),
`glmnet` (penalized GLM benchmark), `mixOmics` (sparse / group PLS).
