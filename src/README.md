# SSDD — Raw Building-Level Metrics

Computes two raw per-building geometric metrics from a footprint layer:
**SD** (structure density) and **SS** (structure separation). Optionally
spatial-joins CAL FIRE DINS damage-inspection points. Stops there on purpose —
normalization, blending and predictive modeling are downstream analysis
choices, not package behavior.

The metric set and functional forms were selected by the metric-specification
experiments in `dev/py/metric-specification/`, which screened the earlier
four-metric set {KD, BA, DP, OP} and functional-form variants on univariate
separation of DINS destroyed/survived labels across the Eaton, Palisades and
Mountain fires. SD replaces KD + BA (the weighting choice subsumes the
basal-area idea); SS replaces DP + OP (nearest-neighbour distance dominates
every averaged or orientation-weighted alternative).

Both metrics are coded so **higher = more vulnerable** (denser surroundings /
closer nearest neighbour).

---

## What the package computes

### `SD` — structure density

Kernel-weighted neighbour mass per unit area within `r_D` of each building's
centroid:

```
SD_i = (1 / (pi * r_D^2)) * sum_{j != i, d_ij <= r_D} w_j * K(d_ij / r_D)
```

- `d_ij` — centroid-to-centroid distance (meters).
- `K` — beta-family kernel `K_n(u) = (1 - u^2)^n`: `uniform` (n=0, default),
  `epanechnikov`, `quartic`, `triweight`. Kernel shape is immaterial to
  separation (all forms within ~0.02 AUC), so the simplest is the default.
- `w_j = sqrt(area_j)` by default (`weight="unit"` → `w_j = 1`, a smoothed
  count; `weight="area"` → `w_j = area_j`, a smoothed basal area).
- Self is excluded; an isolated building scores 0.

### `SS` — structure separation

Inverse wall-to-wall distance to the (true) nearest neighbour:

```
SS_i = 1 / (d_nn + epsilon)          d_nn = min_{j != i} d_ij
```

- `d_ij` — polygon-to-polygon (wall-to-wall) distance, so touching walls give
  `d = 0` and SS saturates at `1 / epsilon`.
- **No isolation cutoff**: `nn` uses the true nearest neighbour and ignores
  `r_S`, so a remote structure keeps its distance ordering (`1/(d+ε)`) instead
  of collapsing to 0. Only a structure alone in the layer scores 0.
- No averaging over neighbours and no orientation weighting **by default**: the
  specification experiments found the plain nearest-neighbour form dominates
  every averaged / orientation-weighted alternative (pooled AUC 0.619 vs ≤ 0.59).
  The alternatives remain reachable via `agg` (`nn` / `uniform` / `power1` /
  `power2`) and `orient` (`flat` / `gauss` / `cos2` / `cos4`) — e.g. the legacy
  DP is `agg="power1"`, the legacy OP is `agg="power1", orient="gauss"`. `r_S`
  bounds the neighbour set for the **averaging aggs only**; `nn` disregards it.

### Attributes carried alongside

| column | meaning |
|---|---|
| `ssdd_id` | stable integer id (added if missing) |
| `bld_area` | footprint area (m²) |
| `phi_deg` | dominant footprint orientation — longest edge of the minimum rotated rectangle, in [0, 180). Not a metric input; kept for downstream analyses. |
| `cent_x`, `cent_y` | true polygon centroid coordinates (analysis CRS) |

---

## Interpretation cheat sheet

| situation | SD | SS |
|---|---|---|
| remote structure, nearest 200 m away | ~0 | 1/200.5 ≈ 0.005 |
| exurban lot, one neighbour 30 m away | small | 1/30.5 ≈ 0.033 |
| suburban street, nearest wall 10 m | moderate | 1/10.5 ≈ 0.095 |
| dense infill, touching row-houses | large | 1/ε = 2.0 (cap) |

SS is a per-structure quantity (its variance is almost entirely within
neighbourhoods); SD is a neighbourhood quantity (half or more of its variance
is between 600 m blocks). They are complementary, not redundant.

---

## Inputs

### Building footprints (required)

Any OGR-readable polygon layer (SHP / GPKG / GeoJSON). Reprojected internally
to the target EPSG (default 32611, UTM 11N) — geometry must end up in meters.

### DINS structure points (optional)

Point layer with damage attributes (e.g. `DAMAGE`). Joined point-in-polygon
onto footprints (`how="left"` keeps all buildings; `--dins-only` keeps only
matched ones).

### Gotcha: identifier column

`join_dins` de-duplicates on `ssdd_id` when multiple DINS points fall in one
footprint (keeps the first). If your layer already has a column named
`ssdd_id` it will be respected, not overwritten.

---

## Tuning knobs

| param | default | affects | notes |
|---|---|---|---|
| `r_D` | 200 m | `SD` | Neighbourhood radius. Separation rises monotonically over the swept grid (pooled AUC 0.646 @50 → 0.706 @200); 200 is the top of the swept range, so larger radii are untested. |
| `r_S` | 50 m | `SS` | Neighbour window for the **averaging** aggs (`uniform`/`power1`/`power2`) only. The default `nn` uses the true nearest neighbour and ignores it. |
| `epsilon` | 0.5 m | `SS` | Distance floor; bounds SS at 1/ε for touching walls. Held fixed in the specification experiments — not empirically tuned. |
| `kernel` | `"uniform"` | `SD` | `uniform` / `epanechnikov` / `quartic` / `triweight`. Shape immaterial; simplest preferred. |
| `weight` | `"root_area"` | `SD` | `unit` (count), `area` (basal area), `root_area`. Best pooled separation; never worse than 2nd in any single fire. |
| `agg` | `"nn"` | `SS` | Neighbour aggregation: `nn` (nearest only), `uniform` (mean g), `power1`/`power2` (mean g/(d+ε)^p). nn is the specified form. |
| `orient` | `"flat"` | `SS` | Orientation weight g(θ): `flat` (g=1), `gauss`, `cos2` (Lambert), `cos4`. flat = no orientation preference. |
| `sigma` | 10° | `SS` | Gaussian orientation tolerance; only used when `orient="gauss"`. |

### Behavioral notes

- **SD truncates hard at `r_D`** — a neighbour at `r_D + 0.1` contributes 0.
- **SS ignores everything but the nearest wall.** Adding a second neighbour
  changes SD but not SS.
- **`epsilon` sets the SS ceiling.** All shared-wall pairs collapse to exactly
  `1/epsilon`; degrees of "touching" are not distinguished.
- **SD is 0 for a building with no neighbour in `r_D`; SS is 0 only for a
  building alone in the layer** (nn has no isolation cutoff). Both zeros sit at
  the *low-risk* end by the design story.

---

## Outputs

### `{run_name}_raw_metrics.csv`

One row per building: `ssdd_id`, `bld_area`, `phi_deg`, `cent_x`, `cent_y`,
`SD`, `SS`, then any DINS columns from the join.

### `{run_name}_buildings.gpkg`

Same table with geometry, layer `buildings_raw`.

### `{run_name}_compute_log.txt`

Run parameters, input paths, counts, elapsed time.

---

## Installation

```bash
conda activate ssdd          # geopandas, shapely>=2, pyproj, tqdm, pandas
cd SSDD/src
```

No build step — the package is imported from `src/` (add it to `sys.path` or
run the CLI from `src/`).

---

## Programmatic API

```python
import sys; sys.path.insert(0, "path/to/SSDD/src")

from ssdd import (
    read_buildings, ensure_projected_meters,
    RawMetricParams, compute_raw_metrics,
)

bld = read_buildings("buildings.gpkg")
bld = ensure_projected_meters(bld, 32611)

params = RawMetricParams(r_D=200.0, r_S=50.0, epsilon=0.5,
                         kernel="uniform", weight="root_area")
out = compute_raw_metrics(bld, params=params)
out[["ssdd_id", "SD", "SS"]]
```

### Per-metric handles

```python
from ssdd import compute_SD_series, compute_SS_series

sd = compute_SD_series(bld, r_D=200.0)                 # Series named SD
ss = compute_SS_series(bld, r_S=50.0, epsilon=0.5)     # Series named SS
```

Both accept prebuilt STRtrees / geometry arrays to avoid rebuilding indexes
when computing several variants.

### Geometry helpers

```python
from ssdd import kernel_value, dominant_orientation_degrees, angle_difference_deg

kernel_value(0.5, "uniform")            # 1.0
dominant_orientation_degrees(poly)      # long-axis azimuth in [0, 180)
angle_difference_deg(10.0, 170.0)       # 20.0 — folded to [0, 90]
```

---

## Command-line interface

```bash
cd SSDD/src
python ssdd_compute.py \
    --buildings ../_data/raw/buildings/LARIAC6_Buildings_2020_eaton.shp \
    --dins      ../_data/raw/dins/DINS_2025_Eaton_Public_View.geojson \
    --output    ../_data/processed/eaton \
    --run-name  eaton
```

- `--r-d`, `--r-s`, `--epsilon`, `--kernel`, `--weight` — override the matching
  `RawMetricParams` field.
- `--dins-only` — inner join (keep only DINS-matched buildings).
- `--epsg` — analysis CRS (default 32611).

---

## Synthetic geometry fixtures

`ssdd.synthetic` builds parametric scenes (isolated building, pairs at
controlled spacing/orientation, grids, heterogeneous-size clusters, L-shapes) in a projected CRS, ready for
`compute_raw_metrics`. Used by the test suite, where every expectation is
hand-computed from the formulas above.

| generator | scene |
|---|---|
| `isolated_building(...)` | one rectangle, no neighbours |
| `pair(spacing, orientation_offset_deg, …)` | two rectangles, controlled wall gap and relative rotation |
| `touching_pair(...)` | shared wall (SS saturation case) |
| `grid(n, pitch, …)` | n×n regular lattice (SD truncation counts) |
| `sized_cluster(specs)` | explicit `(cx, cy, w, l[, angle])` rectangles — heterogeneous size / orientation, hand-checkable |
| `lshape(arm_length, arm_width, angle_deg, …)` | single L-shaped (concave) footprint — non-rectangular MRR / wall-to-wall |
| `random_cloud(n, …, area_range, angle_range)` | non-overlapping cloud; optional size / rotation variation for sensitivity |

---

## Tests

```bash
cd SSDD
PYTHONPATH=src python -m pytest tests/ -q
```

Closed-form checks of SD (kernel family, weight variants, truncation counts)
and SS (inverse distance, saturation at 1/ε, r_S truncation, orientation
invariance) on synthetic scenes.
