# SSDD

This repository holds the project code for the **Structure Separation
Distance Density (SSDD)** work. It's an actively developing project; the
layout below documents what's currently checked in, not what may be added
later.

## What's in the repo right now

- A Python package (`src/ssdd/`) that computes the four raw per-building
  SSDD metrics — KD, BA, DP, OP — from building footprints and, optionally,
  CAL FIRE DINS structure-point inputs.
- A command-line driver (`src/ssdd_compute.py`) and a batch wrapper
  (`scripts/run_la_fires.sh`) for processing fires through the package.
- A test suite (`tests/`) of analytic correctness checks on synthetic
  geometries.
- The original SSDD notebook and supporting reference material under `dev/`.

Algorithmic detail, function signatures, units / ranges, and tuning
behavior for the package live in [`src/README.md`](src/README.md). The rest
of this file describes how the repository is organized and how to run the
LA-fires workflow.

## Repo layout

```
SSDD/
├── src/                  Package source + CLI driver
│   ├── README.md         Algorithm details, tuning knobs, programmatic API
│   ├── ssdd_compute.py   CLI for computing raw metrics on a footprint file
│   └── ssdd/             Importable package (after pip install -e .)
├── tests/                pytest suite — analytic checks on synthetic geometries
├── scripts/              Batch wrappers (e.g. run_la_fires.sh)
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
python -m pytest tests/ -v          # 12 tests, should all pass
```

That installs the `ssdd` package in editable mode. See
[`src/README.md`](src/README.md) for how to call it programmatically or
through the `src/ssdd_compute.py` CLI.

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

- The Eaton and Palisades building-footprint and DINS data used here are the
  same inputs prepared for Kenny, Johns, Pawlak, Fricker, Yost, & Ritter
  (2026), *Urban trees and structure loss in the 2025 Eaton and Palisades
  fires*, **Urban Forestry & Urban Greening** 121,
  [doi:10.1016/j.ufug.2026.129470](https://doi.org/10.1016/j.ufug.2026.129470).
- **LARIAC6 building footprints** (`LARIAC6_Buildings_2020_*.shp`) come from
  the LA Region Imagery Acquisition Consortium; obtain from the dataset
  owner. Either capitalization of the fire-name suffix (`eaton` vs `Eaton`)
  is acceptable on a case-insensitive filesystem.
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
convex-blending steps that the `ssdd` package intentionally omits — lives
in [`dev/notebooks/SSDD.ipynb`](dev/notebooks/SSDD.ipynb). A single-file
Python script derived from that notebook is in
[`dev/scripts/py/ssdd.py`](dev/scripts/py/ssdd.py). Neither is loaded by
the package in `src/`.
