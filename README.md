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
├── dev/                  Historical / exploratory material — not loaded by src
│   ├── notebooks/        Original SSDD prototype notebook
│   └── scripts/py/       Single-file script derived from the notebook
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

For the LA fires (Palisades + Eaton) batch, see
[`scripts/run_la_fires.sh`](scripts/run_la_fires.sh).

## Editor / type-checker setup

`pyrightconfig.json` is configured for the type-checker that Positron and
VS Code's Pylance share. After cloning, install the package once
(`pip install -e .`) and select the `ssdd` conda env as the workspace
interpreter; imports and type checks should resolve cleanly.

## Reference

The convex-blending scheme that originally defined SSDD is preserved as
historical reference in
[`dev/scripts/py/ssdd.py`](dev/scripts/py/ssdd.py) (derived from
[`dev/notebooks/SSDD.ipynb`](dev/notebooks/SSDD.ipynb)). The production
package in `src/` does not depend on it.
