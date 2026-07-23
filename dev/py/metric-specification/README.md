# Metric specification

Justifies the SSDD metric design in two steps: **how many metrics** (collapse the
four raw metrics KD/BA/DP/OP into two, SD/SS) and **which functional form** each
of the two should take. Both are decided on *univariate* separation — each metric
is assessed alone against the DINS destroyed/survived labels, at two spatial
scales (individual structure, and 600 m neighborhood).

Four numbered scripts run in order from the repo root (or cell-by-cell in
Positron — each has `# %%` cells and resolves paths by walking up to `_data/`),
plus a standalone footprint-figure script.

| step | script | reads | writes |
|---|---|---|---|
| 0 | `00_prepare.py` | raw footprints + DINS | `_data/processed/{fire}/{fire}_buildings.gpkg` |
| 1 | `01_primitives.py` | `_data/processed/{fire}/*_buildings.gpkg` | `_out/focal.parquet`, `_out/pairs.parquet` |
| 2 | `02_metrics.py` | primitives | `_out/metrics.parquet` (wide), `_out/forms.csv` |
| 3 | `03_compare.py` | metrics | `_out/compare_primary.csv`, `_out/compare_profiles.csv`, `_out/compare_collapse.csv` |
| 4 | `04_plots.py` | comparison CSVs | `_out/*.png` (pooled figures + per-fire tables) |

Footprint maps (block-grid overlays + bare footprints) are **not** part of this
pipeline — they are figures about the data, not the metrics — and live in
`../plot-footprints.py`, writing to `dev/py/_out/`.

```bash
PY=~/anaconda3/envs/ssdd/bin/python
for s in 01_primitives 02_metrics 03_compare 04_plots; do $PY dev/py/metric-specification/$s.py; done
$PY dev/py/plot-footprints.py       # footprint figures (independent of this pipeline)
```

## The two metrics (design grids in `02_metrics.py`)

**SD — structure density** (collapses KD + BA):
`SD = Σ_j w_j · K(d_cent_j / r_D) / (π r_D²)`, neighbors within `r_D`.
- kernel `K`: uniform · sqrt · Epanechnikov · quartic · triweight · exponential
- weight `w`: `unit` (→ smoothed count = KD) · `area` (→ smoothed basal area = BA) · `root_area`
- `r_D ∈ {50, 100, 150, 200}`

**SS — structure separation** (collapses DP + OP): orientation-weighted proximity
to neighbors within `r_S`.
- orientation `g(θ)`: `flat` (g≡1 = DP) · gaussian σ∈{5,10,20} · cos⁴θ · Lambert cos²θ (= OP)
- aggregation: `nn` (nearest only) · `uniform` (mean g, no distance) · `power1`
  (mean g/(d+ε)) · `power2` (mean g/(d+ε)²)
- `r_S ∈ {25, 50, 75}`

Radius grids are **univariate**: SD is swept over `r_D` only, SS over `r_S` only —
they are never crossed. Each concept is a *branch* of one grid: KD = SD·unit,
BA = SD·area, DP = SS·flat, OP = SS·(g≠flat). So SD's grid ⊇ {KD, BA} and SS's
grid ⊇ {DP, OP} by construction, which is what makes the collapse checkable.

## Two separation measures (`03_compare.py`)

- **struct_auc** — structure-level AUC, `P(metric ranks a destroyed structure
  above a survived one)`. Rank-based, base-rate free.
- **nbhd_pauc** — neighborhood pseudo-AUC: over pairs of 600 m blocks with
  different damage rates, `P(the form's block-mean orders them the same way as
  damage) + ½ ties`. Somers' D mapped to [0,1]; 0.5 = chance, same scale as AUC.
  Reported with the Spearman correlation between block-mean and block damage rate.

Blocks are 1000 m→**600 m** UTM cells; neighborhood measures use blocks with ≥20
buildings. `compare_profiles.csv` additionally splits each form's variance into
between/within-block (`icc_between` = the neighborhood share).

---

## Comparison 1 — Number of metrics (4 → 2)

**Question.** Does replacing {KD, BA, DP, OP} with {SD, SS} lose univariate
separation?

**Design.** For each concept, take the best form achievable in its branch of the
grid (best by `struct_auc` and by `nbhd_pauc` separately), per fire and pooled.
Because SD nests KD/BA and SS nests DP/OP, the check is whether the *combined*
family's best is ≥ each constituent's best:

> **SD.best ≥ max(KD.best, BA.best)** and **SS.best ≥ max(DP.best, OP.best)**,
> on both scales, in each fire.

The real production metrics (`KD_raw`/`BA_raw`/`DP_raw`/`OP_raw`, carried straight
off the gpkgs at their shipped radii, r_D=100 / r_S=50) are shown alongside as a
reference for where the current package sits — but they mix the collapse with
re-tuning, so the collapse claim rests on the best-in-branch bars, not those dots.

**Outputs.** `compare_collapse.csv` (per fire); `collapse.png` (pooled bars +
production dots); `table_num_metrics.png` (per-fire table).

**Caveat.** Univariate separation shows each new metric ≥ its two constituents
*individually*; it does not prove SD captures the *joint* KD+BA information. That
would need a multivariate fit (RF + spatial CV). The nesting argument plus this
dominance table is the intended justification; the RF check is the upgrade path.

## Comparison 2 — Functional forms (per metric)

**Question.** Which kernel/weight (SD) and orientation/aggregation (SS), at which
radius, gives the best univariate separation — and do the fires agree?

**Design.** Score every form in each grid (SD: 6 kernels × 3 weights × 4 r_D;
SS: 6 orientations × 4 aggregations × 3 r_S) on both scales, per fire and pooled.
Read the heatmaps for the pooled winner and the tables for cross-fire agreement.

**Outputs.** `compare_primary.csv` (every form × fire × scale — the full result);
pooled heatmaps `sd_heat_{struct,nbhd}.png`, `ss_heat_{struct,nbhd}.png`;
`scale_scatter.png` (struct vs nbhd, all forms); `profiles_icc.png` (where each
form's variance lives); per-fire tables `table_forms_sd.png`, `table_forms_ss.png`
(shown at the pooled-winning radius; full radius sweep is in the CSV). Kernel
shapes: `sd_kernels.png`, `ss_orient.png`, `ss_distance.png`.

Plots are **pooled**; the per-fire breakdown is in the `table_*` PNGs and the CSVs.

---

## Scaling to more fires

The pipeline is built to grow. To add fires:

1. put their processed gpkgs at `_data/processed/{fire}/{fire}_buildings.gpkg`
   (same columns as the existing ones: `geometry`, `phi_deg`, `cent_x/y`,
   `bld_area`, `DAMAGE`, and `KD_raw`/`BA_raw`/`DP_raw`/`OP_raw`);
2. add the fire name to `FIRES` at the top of `01_primitives.py`;
3. rerun `01 → 02 → 03 → 04` and `plot-footprints.py`.

Everything downstream **discovers the fire set from the data** (03/04/plots read
it off `focal.parquet` / the CSVs), so tables gain columns and per-fire figures
gain panels automatically — no other edits.

## Package independence

**No code dependency on the `ssdd` package, end to end** — the two helpers needed
(`angle_difference_deg`, `dominant_orientation_degrees`) are inlined in
`01_primitives.py` and `00_prepare.py`. `00_prepare.py` builds the processed gpkgs
(area / centroid / orientation / DINS join) itself, so adding a new fire never
touches the package. It survives package refactors and outlives the metric
implementation it is helping to specify.

## Spatial blocking

Neighborhood measures use a **per-fire adaptive** cell size: `03_compare.py` picks,
from `{600…4000} m` (all ≥3× the largest `r_D`), the size giving the most usable
(≥20-building) blocks — fine for dense fires (Eaton/Palisades 600 m), coarse for
sparse ones (Mountain 1000 m). Very small fires stay block-limited regardless, so
read their `nbhd_pauc` with care.

## Notes

- Step 1 is the slow one (polygon distances); its parquet outputs are cached, so
  rerun it only when the building/DINS data changes.
- The old per-metric slide scripts were folded into `04_plots.py`
  (`sd_kernels`/`ss_orient`/`ss_distance` replace the KD/DP/OP kernel plots) and
  `plot-footprints.py` (block maps + bare footprints). Pre-collapse exploratory
  scripts were removed.
