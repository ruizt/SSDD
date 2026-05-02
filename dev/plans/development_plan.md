# SSDD Development Plan

**Project:** Structure Separation Distance Density  
**Context:** Wildland-Urban Interface (WUI) building risk metrics  
**Date:** 2026-05-02

---

## Overview

SSDD quantifies two complementary properties of each building footprint in a WUI dataset:

- **SD (Structure Density):** How densely a building is surrounded, combining a quartic kernel density estimate (KD) and a basal area fraction (BA) within a fixed-radius window.
- **SS (Separation):** How close and directionally aligned neighboring buildings are, combining a mean inverse wall-to-wall distance proxy (DP) and an orientation-weighted version (OP).
- **SSDD:** A final linear (or geometric) blend of SD and SS, normalized to [0, 1].

The current implementation lives in `dev/notebooks/SSDD.ipynb` — a single Jupyter notebook that loads data, runs all spatial computations, normalizes, blends, and exports results. The goal of this plan is to split that into a reproducible, modular pipeline with clearer ownership between the spatial compute stage and the analytical/blending stage.

---

## Design Principles

1. **Separation of concerns.** Spatial computation (geometry-heavy, slow, Python/GeoPandas) is separate from normalization and blending (fast, parameter-sensitive, R-friendly).
2. **The raw metrics are the contract.** Once `KD_raw`, `BA_raw`, `DP_raw`, `OP_raw` are written out per building, all downstream work is pure table arithmetic — no geometry required.
3. **Parameters belong to the analytical stage.** `alpha_D`, `alpha_S`, `beta`, `NORM_METHOD`, `P_LOW/P_HIGH` are blending/normalization knobs. They should be easy to sweep in R without re-running the costly spatial stage.
4. **Geometry is a join, not a dependency.** The spatial stage also writes a geometry sidecar (GeoPackage keyed on `bld_id`) so results can be re-attached for export or visualization at any point.

---

## Proposed Architecture

```
Input
  └── building footprints (.shp / .gpkg)
        │
        ▼
[Stage 1: Spatial Compute — Python]
  ssdd_compute.py
  ├── Load + clean + reproject
  ├── Compute dominant orientation (phi_deg)
  ├── Compute KD_raw  (kernel density, point-based)
  ├── Compute BA_raw  (basal area fraction, polygon-based)
  ├── Compute DP_raw  (mean inverse wall-to-wall distance)
  ├── Compute OP_raw  (orientation-weighted inverse distance)
  └── Write:
        raw_metrics.csv     (bld_id, KD_raw, BA_raw, DP_raw, OP_raw, SS_neighbors)
        buildings_geom.gpkg (bld_id + geometry only, no metrics)
        compute_log.txt     (input path, CRS, n buildings, elapsed time)
        │
        ▼
[Stage 2: Normalize + Blend — R or Python]
  ssdd_blend.R  /  ssdd_blend.py
  ├── Read raw_metrics.csv
  ├── Normalize KD, BA, DP, OP  (configurable method + percentile bounds)
  ├── Blend SD = alpha_D * KD + (1 - alpha_D) * BA
  ├── Blend SS = alpha_S * DP + (1 - alpha_S) * OP
  ├── Blend SSDD = beta * SD + (1 - beta) * SS
  ├── QA plots (distributions, scatter SD vs SS, spatial maps via join)
  └── Write:
        {run_name}_SSDD_{run_tag}.csv       (full attribute table)
        {run_name}_SSDD_{run_tag}.gpkg      (joined back to geometry)
        {run_name}_SSDD_{run_tag}_summary.txt
```

Stage 2 has no geometry dependency — all inputs are a flat CSV and an optional geometry sidecar for
output. This makes it straightforward to implement in either R or Python; both are first-class options.
The choice can be made per-use-case (e.g., R for interactive exploration and parameter sweeps,
Python for scripted batch runs or integration into a larger pipeline).

---

## File Structure (target)

```
SSDD/
├── SSDD.Rproj
├── README.md
├── .gitignore
│
├── dev/                          ← working/exploratory material
│   ├── notebooks/
│   │   └── SSDD.ipynb            ← original prototype notebook
│   ├── plans/
│   │   └── development_plan.md   ← this file
│   └── scripts/                  ← draft scripts purled from notebook
│       ├── py/ssdd.py
│       └── r/ssdd.R
│
└── src/                          ← production pipeline (to be built)
    ├── py/
    │   ├── ssdd_compute.py       ← Stage 1: spatial compute → raw CSV
    │   └── ssdd_blend.py         ← Stage 2 (Python): normalize + blend + export
    └── r/
        └── ssdd_blend.R          ← Stage 2 (R): normalize + blend + export
```

Both Stage 2 implementations are maintained in parallel and should produce identical numerical outputs
given the same inputs and parameters. Either can be used depending on workflow preference.

Data files (building footprints, outputs) are not committed to the repo and should be stored externally (see `.gitignore`).

---

## Stage 1 — `ssdd_compute.py`

### Inputs (CLI or config block at top of file)
| Parameter | Description | Default |
|---|---|---|
| `--input` | Path to building footprints (SHP/GPKG/GeoJSON) | required |
| `--layer` | Layer name for GPKG/FGDB | `None` |
| `--output` | Output directory | required |
| `--run-name` | Prefix for output files | `"run"` |
| `--epsg` | Target CRS (equal-area, meters) | `3310` |
| `--r-d` | SD buffer radius (m) | `100.0` |
| `--r-s` | SS search radius (m) | `50.0` |
| `--epsilon` | Distance floor (m) | `0.5` |
| `--sigma-theta` | Orientation tolerance (deg) | `15.0` |
| `--weight-by-area` | Weight KD neighbors by footprint area | `False` |

### Outputs
- `raw_metrics.csv` — one row per building:  
  `bld_id, KD_raw, BA_raw, DP_raw, OP_raw, SS_neighbors, phi_deg, bld_area`
- `buildings_geom.gpkg` — geometry-only GeoPackage keyed on `bld_id`
- `compute_log.txt` — input path, CRS, n buildings, parameters used, elapsed time

### Notes
- `r_D`, `r_S`, `epsilon`, `sigma_theta` are spatial parameters and belong here — they define the neighborhood and thus affect the raw values.
- `alpha_D`, `alpha_S`, `beta`, and normalization settings do **not** belong here — they affect only blending.
- The script should be fully re-runnable; if raw outputs already exist and inputs/parameters haven't changed, a `--force` flag should control re-computation.

---

## Stage 2 — `ssdd_blend.R` / `ssdd_blend.py`

Stage 2 has no geometry dependency. Both implementations share the same interface and should produce
identical outputs. Use whichever fits the workflow.

### Inputs (config block at top of file)
| Parameter | Description | Default |
|---|---|---|
| `RAW_CSV` | Path to `raw_metrics.csv` from Stage 1 | required |
| `GEOM_GPKG` | Path to `buildings_geom.gpkg` from Stage 1 | required |
| `OUT_DIR` | Output directory | required |
| `RUN_NAME` | Prefix for output files | `"run"` |
| `alpha_D` | KD/BA blend weight | `0.5` |
| `alpha_S` | DP/OP blend weight | `0.5` |
| `beta` | SD/SS blend weight | `0.5` |
| `NORM_METHOD` | `"minmax"` or `"robust"` | `"robust"` |
| `P_LOW`, `P_HIGH` | Percentile bounds for robust normalization | `2, 98` |

### Outputs
- `{run_name}_SSDD_{run_tag}.csv` — full attribute table (raw + normalized + blended)
- `{run_name}_SSDD_{run_tag}.gpkg` — geometry + all attributes joined
- `{run_name}_SSDD_{run_tag}_summary.txt` — parameters and summary stats

### Notes
- Because blending is cheap, this script is designed to be run repeatedly with different parameter combinations — a parameter sweep function should be straightforward to add.
- QA plots (distributions, SD vs SS scatter, spatial maps) live here.
- The geometric mean variant (`SSDD_geom`) is also computed here.
- Geometry rejoining (`buildings_geom.gpkg` + CSV → GeoPackage) requires `sf` in R or `geopandas` in Python, but is a simple join — no spatial computation involved.

---

## Development Phases

### Phase 1 — Refactor Python (spatial stage)
- [ ] Write `src/py/ssdd_compute.py` from `dev/scripts/py/ssdd.py`
- [ ] Strip out blending/normalization parameters (alpha, beta, NORM_METHOD)
- [ ] Output: `raw_metrics.csv` + `buildings_geom.gpkg` + `compute_log.txt`
- [ ] Add `--force` flag; print clear progress + timing
- [ ] Test on SLO sample dataset

### Phase 2 — Build analytical stage (normalize + blend)
- [ ] Write `src/r/ssdd_blend.R` from `dev/scripts/r/ssdd.R`
- [ ] Write `src/py/ssdd_blend.py` from `dev/scripts/py/ssdd.py` (blend/norm section only)
- [ ] Both: read raw CSV, implement norm + blend, join to geometry for spatial output
- [ ] Both: add QA plots (distributions + spatial map)
- [ ] Cross-validate: R and Python outputs match to numerical tolerance with same inputs + parameters
- [ ] Test round-trip: outputs match notebook results with same parameters

### Phase 3 — Validation
- [ ] Run both stages end-to-end on SLO sample; confirm SSDD output matches notebook
- [ ] Document any numerical differences (e.g., due to `st_distance` vs Shapely `distance`)
- [ ] Update README with usage instructions

### Phase 4 — Extensions (future / optional)
- [ ] Parameter sweep helper in R (grid over alpha_D, alpha_S, beta)
- [ ] Parallel execution in Python (multiprocessing for KD/BA/SS loops)
- [ ] Support for multiple input datasets / batch mode
- [ ] Quarto report template for run summaries

---

## Data Requirements

The following files are **not** in the repo and must be obtained externally:

| File | Description | Where to get |
|---|---|---|
| `Buildings.shp` (or `.gpkg`) | Building footprints for analysis area | Request from repo owner |

The sample used during notebook development covers San Luis Obispo County (~22,974 buildings), originally in EPSG:2874 (California State Plane Zone 5), reprojected to EPSG:3310 for analysis.

---

## Dependencies

**Python — Stage 1 (spatial compute)**
- `geopandas`, `shapely >= 2.0`, `numpy`, `pandas`, `tqdm`
- Conda environment name in notebook: `ssdd`

**Python — Stage 2 (blend, optional)**
- `pandas`, `numpy`, `matplotlib`
- `geopandas` only if writing spatial output

**R — Stage 2 (blend, optional)**
- `sf >= 1.0` (only needed if writing spatial output; `st_minimum_rotated_rectangle` is Stage 1 only)
- `ggplot2`, `pbapply`
