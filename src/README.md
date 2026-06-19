# SSDD — Raw Building-Level Metrics

This package computes four raw, un-blended per-building metrics from a
footprint layer. Everything downstream — normalization, blending, predictive
modeling — is **out of scope**.

The original SSDD definition (including the normalization and convex-blending
steps not implemented here) is preserved in
[`dev/notebooks/SSDD.ipynb`](../dev/notebooks/SSDD.ipynb), with a derived
single-file Python script at
[`dev/scripts/py/ssdd.py`](../dev/scripts/py/ssdd.py). Neither is loaded by
this package.

---

## What the package computes

Four raw metrics per building, in two conceptual families:

### Structure Density (SD)

> "How much building is around me?"

| Metric | What it counts | Distance used | Includes self? |
|---|---|---|---|
| **`KD_raw`** | Kernel-weighted neighbor count per m² | Centroid-to-centroid | No |
| **`BA_raw`** | Fraction of nearby ground covered by buildings | Polygon ∩ buffer overlap | Yes |

### Structure Separation (SS)

> "How close (and how aligned) are my neighbors?"

| Metric | What it counts | Distance used | Filters |
|---|---|---|---|
| **`DP_raw`** | Mean `1/(d + ε)` over neighbors | Wall-to-wall | None |
| **`OP_raw`** | Mean `g(Δφ) · 1/(d + ε)` over neighbors | Wall-to-wall | Down-weights neighbors not aligned with focal |

The SD family treats buildings as a **continuous field** (counts/area or
covered-area-fraction). The SS family treats buildings as a **set of
neighbors** (per-neighbor inverse-distance terms, averaged).

---

## Formulas

### `KD_raw` — quartic kernel density

$$
\\mathrm{KD}_i = \\frac{1}{\\pi r_D^2} \\sum_{\\substack{j \\ne i \\\\ d_{ij} \\le r_D}} w_j \\cdot K\\!\\left(\\frac{d_{ij}}{r_D}\\right),
\\qquad K(u) = (1-u^2)^2 \\text{ for } u\\in[0,1]
$$

- `d_ij` = distance between rep points (centroids for convex polygons)
- `w_j = 1` by default (`weight_by_area=True` → `w_j = area_j`)
- Self excluded — an isolated building has `KD_raw = 0`

### `BA_raw` — basal area fraction

$$
\\mathrm{BA}_i = \\frac{\\sum_j \\mathrm{area}\\!\\left(P_j \\cap \\mathrm{buffer}(P_i, r_D)\\right)}{\\mathrm{area}\\!\\left(\\mathrm{buffer}(P_i, r_D)\\right)}
$$

- The "window" is `buffer(P_i, r_D)` — a Minkowski sum (rounded rectangle for
  rectangular footprints), not a disk centered at a point.
- Self included — isolated building has `BA_raw = self_area / window_area`,
  a small positive number.
- Output is a dimensionless fraction in `[0, 1]`.

### `DP_raw` — mean inverse wall-to-wall distance

$$
\\mathrm{DP}_i = \\frac{1}{m_i} \\sum_{\\substack{j \\in N_i}} \\frac{1}{d_{ij} + \\epsilon}
$$

- `N_i` = neighbors with wall-to-wall distance `d_ij ≤ r_S`
- `m_i = |N_i|`. If `m_i = 0`, `DP_raw = 0`.
- `ε` is a distance floor; caps the per-term contribution at `1/ε` when walls touch.

### `OP_raw` — orientation-weighted variant

$$
\\mathrm{OP}_i = \\frac{1}{m_i} \\sum_{\\substack{j \\in N_i}} \\frac{g(\\Delta\\varphi_{ij})}{d_{ij} + \\epsilon},
\\qquad g(\\theta) = \\exp\\!\\left(-(\\theta/\\sigma_\\theta)^2\\right)
$$

- `Δφ_ij` = angle between the longest edges of buildings *i* and *j*, folded into `[0°, 90°]`.
- `σ_θ` is the orientation tolerance. Smaller → alignment matters more.
- Always `OP_raw ≤ DP_raw`. Equal only when all neighbors are perfectly parallel.

---

## Interpretation cheat sheet

| Setup | `KD` | `BA` | `DP` | `OP` |
|---|---|---|---|---|
| Isolated building | 0 | tiny (self / window) | 0 | 0 |
| Two close parallel buildings | small | small | large | ≈ DP |
| Two close perpendicular buildings | small | small | large | ≈ 0 |
| Trailer park (many small structures) | high | moderate | high | depends on alignment |
| Mall (one giant footprint) | low | very high | low (or 0 if isolated) | low |
| Dense suburb of similar lots | high | high | high | ≈ DP if streets are gridded |

What it means physically (rough mappings to WUI ignition/spread thinking):

- **`KD`** answers *how many buildings* nearby. Picks up neighborhood character.
- **`BA`** answers *how much built mass* nearby. Picks up the "fuel" density.
- **`DP`** answers *how close the nearest neighbors are*. Inverse-weighted, so the nearest dominate the average.
- **`OP`** answers *how close and aligned* — captures the geometric setup where two parallel walls face each other (relevant for radiative transfer between structures).

---

## Inputs

### Building footprints (required)

- Any GeoPandas-readable layer: Shapefile, GeoPackage, GeoJSON, FGDB, FlatGeobuf.
- `Polygon` or `MultiPolygon` geometries. Null and non-polygonal features are dropped.
- All polygons are passed through `buffer(0)` to repair self-intersections.
- CRS must be defined. It's reprojected to the target CRS before any math.

### DINS structure points (optional)

- Any GeoPandas-readable point layer.
- CRS must be defined; reprojected to the analysis CRS before joining.
- Joined point-in-polygon (`predicate="intersects"`). All DINS attribute columns
  ride along into the output.

### Gotcha: identifier column

The package adds a stable integer column **`ssdd_id`** (not `bld_id`) to avoid
case-insensitive collisions with LARIAC6's existing `BLD_ID` when writing to
GeoPackage (SQLite folds column names case-insensitively).

---

## Tuning knobs

All distances are in **meters**. The target CRS must be projected, units = meters.

| Knob | Default | Affects | Notes |
|---|---:|---|---|
| `r_D` | 100 m | `KD_raw`, `BA_raw` | Neighborhood radius for the SD family. Larger → smoother density signal. |
| `r_S` | 50 m | `DP_raw`, `OP_raw`, `SS_neighbors` | Search radius for the SS family. Defines who counts as a "neighbor". |
| `epsilon` | 0.5 m | `DP_raw`, `OP_raw` | Distance floor. Caps per-term contribution at `1/ε` for touching walls. Match it to your minimum resolvable gap. |
| `sigma_theta` | 15° | `OP_raw` | Orientation tolerance in `g(Δφ)`. At σ=15°, a 30° offset already attenuates by ~98%. |
| `kernel` | `"quartic"` | `KD_raw` | KD kernel shape. Only quartic implemented; one-function extension in `geometry.py` if you want others. |
| `weight_by_area` | `False` | `KD_raw` | If True, bigger neighbors count more in KD. Does not touch BA. |
| `r_NN` | 200 m | `dist_to_nearest_building`, `bearing_to_nearest_building` | Search radius for the nearest-neighbor proximity. Returns `NaN` for both columns when no neighbor lies within `r_NN`. |

### Behavioral notes

- **`r_D` doesn't change KD scale much in uniform density.** Doubling `r_D`
  adds neighbors (roughly ∝ `r_D²`) *and* shrinks the denominator `π r_D²` by
  the same factor. The effect is on **locality vs smoothing**, not absolute
  magnitude.
- **`r_S` does change DP/OP magnitude.** A larger `r_S` admits more distant
  (lower-weighted) neighbors and dilutes the average — DP and OP fall.
- **`ε` saturates inverse-distance terms.** With default `ε = 0.5`, DP and OP
  both top out at 2.0 per term. Smaller `ε` → larger saturation cap →
  touching buildings dominate even more.

---

## Outputs

Three files per run, all written under `--output/`:

### `{run_name}_raw_metrics.csv`

One row per building. Always includes:

| Column | Units | Range | Family | Meaning |
|---|---|---|---|---|
| `ssdd_id` | — | `0 … N-1` | — | Package-assigned integer ID. |
| `bld_area` | m² | `> 0` | — | Footprint area in the analysis CRS. |
| `phi_deg` | ° | `[0, 180)` | — | Angle of the longest MRR edge. |
| `cent_x`, `cent_y` | m | — | — | Representative-point coords in the analysis CRS. |
| `KD_raw` | 1/m² | `≥ 0` | **SD** | Quartic-kernel structure density. |
| `BA_raw` | — | `[0, 1]` | **SD** | Basal area fraction in window. |
| `DP_raw` | 1/m | `[0, 1/ε]` | **SS** | Mean inverse wall-to-wall distance to neighbors. |
| `OP_raw` | 1/m | `[0, 1/ε]` | **SS** | Orientation-weighted `DP_raw`. |
| `SS_neighbors` | — | `≥ 0` | **SS** | Neighbor count within `r_S`. |
| `dist_to_nearest_building` | m | `[0, r_NN]` or NaN | — | Wall-to-wall distance to the nearest other building. NaN if none within `r_NN`. |
| `bearing_to_nearest_building` | ° | `[0, 360)` or NaN | — | Compass bearing (0 = N, clockwise) from focal centroid to nearest-neighbor centroid. NaN if no neighbor within `r_NN`. |

Plus every column from the source footprint layer and every column from the
DINS layer (when `--dins` is supplied). Rows without a DINS hit have NaN in
the DINS columns when `--dins-only` is not set.

### `{run_name}_buildings.gpkg`

Same data + geometry, layer `buildings_raw`. Useful for QGIS inspection, for
spatial joins downstream, or for rejoining metric outputs to footprints after
you've done attribute-only analysis.

### `{run_name}_compute_log.txt`

Plain text: inputs, output paths, parameters used, building count, elapsed.
Enough to reproduce or compare runs.

---

## Installation

Once per environment, from the repo root:

```bash
conda activate ssdd
pip install -e .
```

This installs the package in editable mode — `from ssdd import …` works from
any directory and your local edits to `src/ssdd/` take effect immediately,
without re-installing.

---

## Programmatic API

After installation, from anywhere:

```python
from ssdd import (
    read_buildings, read_dins,
    ensure_projected_meters, join_dins,
    compute_raw_metrics,
)
from ssdd.pipeline import RawMetricParams

bld = read_buildings("_data/raw/buildings/LARIAC6_Buildings_2020_eaton.shp")
bld = ensure_projected_meters(bld, target_epsg=32611)

params = RawMetricParams(r_D=100.0, r_S=50.0, epsilon=0.5, sigma_theta=15.0)
bld = compute_raw_metrics(bld, params=params)

dins = read_dins("_data/raw/dins/DINS_2025_Eaton_Public_View.geojson")
dins = ensure_projected_meters(dins, target_epsg=32611)
bld = join_dins(bld, dins, how="left")   # how="inner" for the DINS-only subset
```

### Per-metric handles

If you want one metric at a time (e.g., for sensitivity sweeps over a single
knob), the metric functions are public:

```python
from ssdd.metrics import (
    compute_KD_series,        # -> pd.Series indexed by gdf.index
    compute_BA_series,
    compute_SS_terms_df,      # -> DataFrame with DP_raw, OP_raw, SS_neighbors
)
```

Each builds its own STRtree if not supplied. Pass `tree_polys=`/`tree_pts=` to
share an index across calls.

### Geometry helpers

```python
from ssdd.geometry import (
    dominant_orientation_degrees,
    angle_difference_deg,
    orientation_factor,
    quartic_kernel,
)
```

Pure functions, no I/O.

---

## Command-line interface

`ssdd_compute.py` is a thin wrapper around `compute_raw_metrics`. After
`pip install -e .` you can run it from anywhere; here, invoked from the repo
root:

```bash
python src/ssdd_compute.py \
  --buildings _data/raw/buildings/LARIAC6_Buildings_2020_eaton.shp \
  --dins      _data/raw/dins/DINS_2025_Eaton_Public_View.geojson \
  --output    _data/processed/eaton \
  --run-name  eaton
```

For the LA-fires batch, see [`scripts/run_la_fires.sh`](../scripts/run_la_fires.sh).

Useful flags:

- `--dins-only` — inner-join with DINS; drops buildings with no DINS hit.
- `--epsg 32611` — change target CRS (default UTM 11N, meters).
- `--r-d`, `--r-s`, `--epsilon`, `--sigma-theta`, `--weight-by-area` — override
  the matching `RawMetricParams` field.

`python ssdd_compute.py --help` for the full list.

---

## Synthetic geometry fixtures

`ssdd.synthetic` provides parametric generators so you can verify the
implementation against hand-computed expected values and run sensitivity
studies on the raw metrics.

| Function | Builds | Validates |
|---|---|---|
| `isolated_building` | One rectangle. | KD/DP/OP = 0; BA = self/window. |
| `pair(spacing, orientation_offset_deg, …)` | Two rectangles at known wall-to-wall distance and relative angle. | Closed-form DP/OP for arbitrary spacing and orientation. |
| `touching_pair` | Two rectangles sharing a wall. | Saturation: `DP_raw = 1/ε`. |
| `grid(n, pitch, …)` | n×n regular lattice. | Neighbor-count partitions and orientation-uniform density. |
| `random_cloud(n, extent, seed)` | n non-overlapping random rectangles. | Benchmarking and realistic-density sweeps. |

All return `GeoDataFrame` in EPSG:32611. Example — sensitivity to `r_D`:

```python
from ssdd import synthetic
from ssdd.pipeline import RawMetricParams, compute_raw_metrics

gdf = synthetic.grid(n=5, pitch=25.0, width=10.0, length=10.0)
for r in (50, 100, 200, 400):
    out = compute_raw_metrics(gdf, params=RawMetricParams(r_D=r), progress=False)
    print(f"r_D={r:>3}   median KD = {out['KD_raw'].median():.2e}   median BA = {out['BA_raw'].median():.3f}")
```

---

## Tests

```bash
conda activate ssdd
cd SSDD
python -m pytest tests/ -v
```

Twelve cases covering:

- **Isolated** — KD/DP/OP = 0; BA = self/window.
- **Parallel pair** — DP = 1/(spacing + ε), OP = DP, KD = `K(d_c/r_D) / (π r_D²)`.
- **Perpendicular pair** — OP attenuated to ~0 by `g(90°)`.
- **Touching pair** — DP and OP saturate at `1/ε`.
- **3×3 grid** — neighbor-count partition (3 corner, 5 edge, 8 center).
- **Orientation sweep** — `OP/DP = g(Δφ)` for Δφ ∈ {0, 15, 30, 45, 60, 75, 90}°.

