"""Raw SSDD metrics: SD (structure density) and SS (structure separation).

These are the two metrics selected by the metric-specification experiments
(``dev/py/metric-specification/``), which screened the earlier four-metric set
{KD, BA, DP, OP} and their functional-form variants on univariate separation
of DINS destroyed/survived labels. SD replaces KD + BA (weighting subsumes the
basal-area idea); SS replaces DP + OP (nearest-neighbour distance dominates
every averaged or orientation-weighted alternative).

Each metric is a standalone function operating on a GeoDataFrame and returning
a pandas Series aligned to ``gdf.index``. Spatial indexes are built internally
if not supplied; pass them in to avoid rebuilding when running both metrics.

Geometry must already be in a projected CRS with units of meters — call
:func:`ssdd.io.ensure_projected_meters` upstream.
"""

from __future__ import annotations

import math
from typing import Any, Optional

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.strtree import STRtree
from tqdm.auto import tqdm

from .geometry import angle_difference_deg, kernel_value

WEIGHTS = ("unit", "area", "root_area")
AGGS = ("nn", "uniform", "power1", "power2")
ORIENTS = ("flat", "gauss", "cos2", "cos4")

# Several parameters here accept "any array-like sequence of geometries / floats"
# — GeometryArray, numpy.ndarray, plain lists. We type them as Any rather than
# enumerate a union, both for readability and to sidestep stub friction between
# geopandas, shapely, and numpy.


def _tree_query_indices(tree: STRtree, query_geom) -> np.ndarray:
    """Return integer indices from STRtree.query(), Shapely 2.x compatible."""
    res = tree.query(query_geom)
    if isinstance(res, (list, tuple, np.ndarray)) and len(res) > 0:
        if isinstance(res[0], (int, np.integer)):
            return np.asarray(res, dtype=int)
    idx_map = {id(g): i for i, g in enumerate(tree.geometries)}
    return np.asarray([idx_map[id(g)] for g in res], dtype=int)


def compute_SD_series(
    buildings: gpd.GeoDataFrame,
    r_D: float,
    kernel: str = "uniform",
    weight: str = "root_area",
    tree_cents: Optional[STRtree] = None,
    cents: Any = None,
    areas: Any = None,
    progress: bool = True,
) -> pd.Series:
    """Structure density: kernel-weighted neighbour mass per unit area.

    SD_i = (1 / (pi * r_D^2)) * sum_{j != i, d_ij <= r_D} w_j * K(d_ij / r_D)

    where d_ij is the centroid-to-centroid distance. ``weight`` selects w_j:

    ``unit``       w_j = 1             — smoothed structure count
    ``area``       w_j = area_j        — smoothed basal area (the old BA idea)
    ``root_area``  w_j = sqrt(area_j)  — default; best pooled separation and
                                          never worse than 2nd in any single fire

    Higher SD = denser surroundings = (design story) more vulnerable.
    """
    if weight not in WEIGHTS:
        raise ValueError(f"Unknown weight {weight!r}; expected one of {WEIGHTS}")
    n = len(buildings)
    if cents is None:
        cents = buildings.geometry.centroid.values
    if tree_cents is None:
        tree_cents = STRtree(cents)
    if weight != "unit" and areas is None:
        areas = buildings.geometry.area.to_numpy()

    out = np.zeros(n, dtype=float)
    norm = math.pi * r_D * r_D
    iterator = range(n)
    if progress:
        iterator = tqdm(iterator, desc="  SD (structure density)")
    for i in iterator:
        ci = cents[i]
        idxs = _tree_query_indices(tree_cents, ci.buffer(r_D))
        total = 0.0
        for j in idxs:
            if j == i:
                continue
            dist = ci.distance(cents[j])
            if dist > r_D:
                continue
            if weight == "unit":
                w = 1.0
            elif weight == "area":
                w = float(areas[j])
            else:  # root_area
                w = math.sqrt(float(areas[j]))
            total += w * kernel_value(dist / r_D, kernel=kernel)
        out[i] = total / norm

    return pd.Series(out, index=buildings.index, name="SD")


def _orient_weight(theta_deg: float, orient: str, sigma: float) -> float:
    """Orientation weight g(theta) for a folded angle difference in [0, 90]."""
    if orient == "flat":
        return 1.0
    if orient == "gauss":
        return math.exp(-((theta_deg / sigma) ** 2))
    if orient == "cos2":
        return math.cos(math.radians(theta_deg)) ** 2
    if orient == "cos4":
        return math.cos(math.radians(theta_deg)) ** 4
    raise ValueError(f"Unknown orient {orient!r}; expected one of {ORIENTS}")


def compute_SS_series(
    buildings: gpd.GeoDataFrame,
    r_S: float,
    epsilon: float,
    agg: str = "nn",
    orient: str = "flat",
    sigma: float = 10.0,
    phi_deg: Any = None,
    tree_polys: Optional[STRtree] = None,
    polys: Any = None,
    progress: bool = True,
) -> pd.Series:
    """Structure separation over neighbours within ``r_S`` (wall-to-wall).

    The default — and the form selected by the metric-specification
    experiments — is the plain nearest-neighbour inverse distance:

    SS_i = g(theta*) / (d* + epsilon),   * = the true nearest neighbour

    ``agg`` selects the aggregation over the neighbour set N_i:

    ``nn``       g(theta*) / (d* + epsilon) at the nearest neighbour (default);
                 uses the TRUE nearest neighbour and ignores ``r_S`` — there is
                 no isolation cutoff, so a remote structure keeps its distance
                 ordering (1/(d+eps)) instead of collapsing to 0.
    ``uniform``  mean_j g(theta_j)                       (no distance decay)
    ``power1``   mean_j g(theta_j) / (d_j + epsilon)     (the legacy DP / OP)
    ``power2``   mean_j g(theta_j) / (d_j + epsilon)^2

    ``r_S`` bounds the neighbour set for the averaging aggregations only
    (``uniform`` / ``power1`` / ``power2``); ``nn`` disregards it.

    ``orient`` selects the orientation weight g on the folded angle difference
    theta in [0, 90] between focal and neighbour ``phi_deg``:

    ``flat``   g = 1 (default — orientation adds no univariate separation)
    ``gauss``  g = exp(-(theta / sigma)^2)
    ``cos2``   g = cos^2(theta)   (Lambert)
    ``cos4``   g = cos^4(theta)

    ``phi_deg`` (per-building dominant orientation) is required whenever
    ``orient != "flat"``. A building alone in the layer (nn) or with no
    neighbour within ``r_S`` (averaging aggs) scores 0; touching walls (d = 0)
    saturate at g / epsilon.

    Higher SS = nearer / more aligned neighbours = (design story) more
    vulnerable.
    """
    if agg not in AGGS:
        raise ValueError(f"Unknown agg {agg!r}; expected one of {AGGS}")
    if orient not in ORIENTS:
        raise ValueError(f"Unknown orient {orient!r}; expected one of {ORIENTS}")
    if orient != "flat" and phi_deg is None:
        raise ValueError("phi_deg is required when orient != 'flat'")

    n = len(buildings)
    if polys is None:
        polys = buildings.geometry.values
    if tree_polys is None:
        tree_polys = STRtree(polys)
    phi = None if phi_deg is None else np.asarray(phi_deg, dtype=float)

    out = np.zeros(n, dtype=float)
    iterator = range(n)
    if progress:
        iterator = tqdm(iterator, desc="  SS (structure separation)")
    for i in iterator:
        Pi = polys[i]

        if agg == "nn":
            # true nearest neighbour (excludes self), no r_S truncation
            nidx, ndist = tree_polys.query_nearest(
                Pi, exclusive=True, all_matches=False, return_distance=True)
            if len(nidx) == 0:
                continue  # alone in the layer -> SS = 0
            j = int(nidx[0])
            d = float(ndist[0])
            g = 1.0 if phi is None else _orient_weight(
                angle_difference_deg(phi[i], phi[j]), orient, sigma)
            out[i] = g / (d + epsilon)
            continue

        idxs = _tree_query_indices(tree_polys, Pi.buffer(r_S))
        total = 0.0
        m = 0
        for j in idxs:
            if j == i:
                continue
            dij = Pi.distance(polys[j])
            if dij > r_S:
                continue
            g = 1.0 if phi is None else _orient_weight(
                angle_difference_deg(phi[i], phi[j]), orient, sigma)
            if agg == "uniform":
                total += g
            elif agg == "power1":
                total += g / (dij + epsilon)
            else:  # power2
                total += g / (dij + epsilon) ** 2
            m += 1
        if m > 0:
            out[i] = total / m

    return pd.Series(out, index=buildings.index, name="SS")
