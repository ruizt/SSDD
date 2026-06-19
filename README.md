# SSDD

**Structure Separation Distance Density** — a Python package that computes
four raw, per-building geometric metrics from a footprint layer, intended as
inputs to downstream Wildland-Urban Interface (WUI) structure-loss modeling.

The package is **strictly scoped to the raw metrics**. Normalization,
blending, and any predictive modeling are left to the user.

## What the package computes

| Family | Metric | Quantifies |
|---|---|---|
| Structure Density (SD) | `KD_raw` | Smoothed neighbor count per m² (kernel density at the focal centroid). |
| Structure Density (SD) | `BA_raw` | Fraction of nearby ground covered by buildings (basal area in a polygon-buffer window). |
| Structure Separation (SS) | `DP_raw` | Mean inverse wall-to-wall distance to neighbors within `r_S`. |
| Structure Separation (SS) | `OP_raw` | `DP_raw` down-weighted for neighbors not aligned with the focal building. |

Full algorithmic and interpretive detail lives in [`src/README.md`](src/README.md).

## Repo layout

```
SSDD/
├── src/                  Package source + CLI driver
│   ├── README.md         Algorithm details, tuning knobs, programmatic API
│   ├── ssdd_compute.py   CLI for computing raw metrics on a footprint file
│   └── ssdd/             Importable package (after pip install -e .)
├── tests/                pytest suite — analytic checks on synthetic geometries
├── scripts/              Reproducibility wrappers (e.g. run_la_fires.sh)
├── dev/                  Reference material — not loaded by src
│   ├── notebooks/        Source notebooks (SSDD.ipynb, overture_wui_buildings.ipynb)
│   └── scripts/py/       Single-file Python script derived from SSDD.ipynb
├── _data/                Local-only inputs and processed outputs (gitignored)
├── pyproject.toml        Build + dependency metadata
├── pyrightconfig.json    Type-checker config for editor integration
└── .vscode/              Positron / VS Code workspace settings
```

## Quick start

```bash
conda env create -n ssdd -c conda-forge python=3.11 \
  geopandas shapely pandas numpy tqdm pytest
conda activate ssdd

cd SSDD
pip install -e .
python -m pytest tests/ -v          # 12 tests
```

To compute raw metrics on a footprint layer, with an optional DINS join:

```bash
python src/ssdd_compute.py \
  --buildings _data/raw/buildings/LARIAC6_Buildings_2020_eaton.shp \
  --dins      _data/raw/dins/DINS_2025_Eaton_Public_View.geojson \
  --output    _data/processed/eaton \
  --run-name  eaton
```

For the LA fires (Palisades + Eaton) batch, see *Processing the LA fires* below.

## Processing the LA fires (Palisades + Eaton)

[`scripts/run_la_fires.sh`](scripts/run_la_fires.sh) is a thin wrapper that
runs `src/ssdd_compute.py` once per fire with consistent file paths and
parameters, producing the raw-metrics + DINS-joined outputs in
`_data/processed/`.

### Required inputs

Stage these files under `_data/raw/` before running. The directory is
gitignored, so files stay local to your machine.

```
_data/raw/
├── buildings/
│   ├── LARIAC6_Buildings_2020_eaton.shp     (+ .shx .dbf .prj sidecars)
│   └── LARIAC6_Buildings_2020_palisades.shp (+ .shx .dbf .prj .cpg sidecars)
└── dins/
    ├── DINS_2025_Eaton_Public_View.geojson
    └── DINS_2025_Palisades_Public_View.geojson
```

Source notes:

- **LARIAC6 building footprints** (`LARIAC6_Buildings_2020_*.shp`) come from
  the LA Region Imagery Acquisition Consortium. The files used here are the
  same ones consumed by Kenny et al. via the LA fires structure-loss
  workflow; obtain from the dataset owner. Either capitalization of the
  fire name suffix (`eaton` vs `Eaton`) is acceptable on a
  case-insensitive filesystem.
- **DINS public-view points** (`DINS_2025_*_Public_View.geojson`) come from
  the CAL FIRE Damage Inspection (DINS) program and are publicly available.

The package reprojects both layers to EPSG:32611 (UTM 11N, meters) at
read time, so input CRS doesn't matter as long as it's defined.

### Running the batch

```bash
conda activate ssdd

# Both fires with default parameters (~50 s for Palisades, ~3 min for Eaton):
./scripts/run_la_fires.sh

# One fire only:
./scripts/run_la_fires.sh palisades

# Override compute parameters (extra flags pass through to ssdd_compute.py):
./scripts/run_la_fires.sh eaton --r-d 150 --r-s 75

# Restrict to the "burned subset" only (inner DINS join):
./scripts/run_la_fires.sh --dins-only
```

The script is location-agnostic — it resolves all paths relative to its
own location, so it works the same whether invoked from the repo root or
from inside `scripts/`. It uses `set -euo pipefail`, so if Palisades fails
Eaton won't run.

### Outputs

For each fire, three files are written under `_data/processed/{fire}/`:

| File | Contents |
|---|---|
| `{fire}_raw_metrics.csv` | One row per building, all attributes incl. `KD_raw`, `BA_raw`, `DP_raw`, `OP_raw`, `phi_deg`, `SS_neighbors`, source-layer columns, and DINS attributes (`DAMAGE`, `STRUCTURETYPE`, …) where they joined. |
| `{fire}_buildings.gpkg` | Same data with footprint geometry, layer `buildings_raw`. Useful for QGIS inspection or downstream spatial joins. |
| `{fire}_compute_log.txt` | Plain-text record of inputs, parameters, building count, and elapsed time — enough to reproduce the run. |

Column definitions, units, ranges and tuning-knob behavior are documented
in [`src/README.md`](src/README.md).

## Editor / type-checker setup

`pyrightconfig.json` is configured for the type-checker that Positron and
VS Code's Pylance share. After cloning, install the package once
(`pip install -e .`) and select the `ssdd` conda env as the workspace
interpreter; imports and type checks should resolve cleanly.

## Reference

The original SSDD definition — including the normalization and
convex-blending steps that this package intentionally omits — lives in
[`dev/notebooks/SSDD.ipynb`](dev/notebooks/SSDD.ipynb). A single-file
Python script derived from that notebook is in
[`dev/scripts/py/ssdd.py`](dev/scripts/py/ssdd.py). Neither is loaded by
the package in `src/`.
