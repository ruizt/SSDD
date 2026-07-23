# dev/r — SSDD tuning & downstream benchmark

Post-processing and evaluation of the Tide `sweep_cv` results. Two scripts, run
in order from the repo root. All outputs land in the gitignored `dev/r/_out/`.

| Script | Does | Key output |
|---|---|---|
| `01_cv_tuning.R` | Treats the CV sweep as a hyperparameter grid: summarises held-out skill across metric form × radius, and selects the tuned config for the **package-default** form and the **optimal** form (each at its best `r_D`), at both scales. | `_out/tuning_selection.csv` + 4 figures |
| `02_benchmark.R` | Ships those tuned configs into a head-to-head against the existing no-SSDD covariates: three structure RFs (`kenny_ext`, `ssdd_default`, `ssdd_optimal`) on identical spatial-block folds. | `_out/benchmark_{results,predictions}.csv` + 2 figures |

```bash
Rscript dev/r/01_cv_tuning.R      # reads _data/processed/sweep_cv/sweep_summary.csv
Rscript dev/r/02_benchmark.R      # reads tuning_selection.csv + sweep_all.csv + Kenny covariates
```

## What the current results say

- **Forms are predictively interchangeable.** Across all 12 form-settings the
  pooled held-out score spans ≤0.02 (structure) / ≤0.04 (neighbourhood) — well
  inside the ~0.07 fold-to-fold noise. The shipped default
  (`SD_uniform_root_area | SS_flat`) sits mid-pack; the "optimal" form beats it
  by ~0.004, i.e. noise.
- **Radius barely matters** over `r_D` 50–300.
- **Fire heterogeneity dominates**: Palisades ≈0.71, Mountain ≈0.58, Eaton ≈0.55
  (structure). Per-fire fits don't rescue Eaton, so its weakness is intrinsic.
- **Tuned SSDD ties extended-Kenny** (Eaton + Palisades, the only fires with
  Kenny building covariates): `ssdd_optimal − kenny_ext ≈ −0.002`,
  `ssdd_default − kenny_ext ≈ −0.009` AUC. SSDD re-encodes the same signal Kenny
  already had — `SS_flat` is the monotone inverse of
  `distance_to_nearest_building` (literal metres), and `SD` parallels
  `build_dens` — so it matches rather than beats those covariates.

## Constraints / notes

- The Kenny head-to-head is **Eaton + Palisades only**. Mountain's covariate
  file is terrain-only (derived from a DEM), with no `build_dens` /
  `area_ext_build` / `distance_to_nearest_building`.
- `02_benchmark.R` mirrors `scripts/tide/sweep_cv/compute.R` exactly (same
  label, folds, ranger settings), so the `ssdd_*` models reproduce the sweep_cv
  per-fire AUCs on the identical sample (eaton=8682, palisades=6041) — up to the
  fact that this fits on 2 fires, not 3.
- The prior 4-metric (KD/BA/DP/OP) RF workspace is archived in the gitignored
  `dev/r/_archive/` (also in git history at commit `6ba0ed6`). The extended-Kenny
  benchmark logic here is ported from that workspace's `benchmark_kenny_ext_rf.R`.
