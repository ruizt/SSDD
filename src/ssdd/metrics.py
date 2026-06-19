"""Raw SSDD metrics: KD, BA, DP, OP.

Each metric is a standalone function operating on a GeoDataFrame and returning
a pandas Series (or DataFrame) aligned to ``gdf.index``. Spatial indexes are
built internally if not supplied; pass them in to avoid rebuilding when
running multiple metrics.

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

from .geometry import (
    angle_difference_deg,
    kernel_value,
    orientation_factor,
)

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


def _representative_points(gdf: gpd.GeoDataFrame) -> Any:
    return gdf.geometry.representative_point().values


def compute_KD_series(
    buildings: gpd.GeoDataFrame,
    r_D: float,
    kernel: str = "quartic",
    weight_by_area: bool = False,
    tree_pts: Optional[STRtree] = None,
    rep_pts: Any = None,
    areas: Any = None,
    progress: bool = True,
) -> pd.Series:
    """Kernel-density-style structure density at each building's rep point.

    KD_i = (1 / (pi * r_D^2)) * sum_{j != i, d_ij <= r_D} w_j * K(d_ij / r_D)
    """
    n = len(buildings)
    if rep_pts is None:
        rep_pts = _representative_points(buildings)
    if tree_pts is None:
        tree_pts = STRtree(rep_pts)
    if weight_by_area and areas is None:
        areas = buildings.geometry.area.to_numpy()

    out = np.zeros(n, dtype=float)
    norm = math.pi * r_D * r_D
    iterator = range(n)
    if progress:
        iterator = tqdm(iterator, desc="  KD (kernel density)")
    for i in iterator:
        ci = rep_pts[i]
        idxs = _tree_query_indices(tree_pts, ci.buffer(r_D))
        total = 0.0
        for j in idxs:
            if j == i:
                continue
            dist = ci.distance(rep_pts[j])
            if dist > r_D:
                continue
            u = dist / r_D
            w = float(areas[j]) if weight_by_area else 1.0
            total += w * kernel_value(u, kernel=kernel)
        out[i] = total / norm

    return pd.Series(out, index=buildings.index, name="KD_raw")


def compute_BA_series(
    buildings: gpd.GeoDataFrame,
    r_D: float,
    tree_polys: Optional[STRtree] = None,
    polys: Any = None,
    progress: bool = True,
) -> pd.Series:
    """Basal-area fraction within a buffer of each footprint.

    BA_i = sum_j area(P_j ∩ buffer(P_i, r_D)) / area(buffer(P_i, r_D))
    """
    n = len(buildings)
    if polys is None:
        polys = buildings.geometry.values
    if tree_polys is None:
        tree_polys = STRtree(polys)

    out = np.zeros(n, dtype=float)
    iterator = range(n)
    if progress:
        iterator = tqdm(iterator, desc="  BA (basal area fraction)")
    for i in iterator:
        win = polys[i].buffer(r_D)
        win_area = win.area
        if win_area <= 0:
            continue
        idxs = _tree_query_indices(tree_polys, win)
        inter_area_sum = 0.0
        for j in idxs:
            inter = polys[j].intersection(win)
            if not inter.is_empty:
                inter_area_sum += inter.area
        out[i] = inter_area_sum / win_area

    return pd.Series(out, index=buildings.index, name="BA_raw")


def compute_SS_terms_df(
    buildings: gpd.GeoDataFrame,
    r_S: float,
    epsilon: float,
    sigma_theta: float,
    phi_deg: pd.Series,
    tree_polys: Optional[STRtree] = None,
    polys: Any = None,
    progress: bool = True,
) -> pd.DataFrame:
    """Separation proxies: mean inverse wall-to-wall distance (DP) and
    orientation-weighted variant (OP), plus neighbor count.

    DP_raw_i = mean_{j in N_i} 1 / (d_ij + eps)
    OP_raw_i = mean_{j in N_i} g(theta_ij) / (d_ij + eps)

    where N_i is the set of buildings within ``r_S`` of building i, d_ij is the
    polygon-to-polygon (wall-to-wall) distance, and g is the orientation weight
    from :func:`ssdd.geometry.orientation_factor`.
    """
    n = len(buildings)
    if polys is None:
        polys = buildings.geometry.values
    if tree_polys is None:
        tree_polys = STRtree(polys)
    phi = phi_deg.to_numpy(dtype=float)

    dp_raw = np.zeros(n, dtype=float)
    op_raw = np.zeros(n, dtype=float)
    m_cnt = np.zeros(n, dtype=int)

    iterator = range(n)
    if progress:
        iterator = tqdm(iterator, desc="  SS (distance + orientation)")
    for i in iterator:
        Pi = polys[i]
        phi_i = float(phi[i])
        idxs = _tree_query_indices(tree_polys, Pi.buffer(r_S))

        inv_sum = 0.0
        inv_orient_sum = 0.0
        m = 0
        for j in idxs:
            if j == i:
                continue
            dij = Pi.distance(polys[j])
            if dij > r_S:
                continue
            inv = 1.0 / (dij + epsilon)
            theta = angle_difference_deg(phi_i, float(phi[j]))
            orient = orientation_factor(theta, sigma_theta)
            inv_sum += inv
            inv_orient_sum += orient * inv
            m += 1

        if m > 0:
            dp_raw[i] = inv_sum / m
            op_raw[i] = inv_orient_sum / m
            m_cnt[i] = m

    return pd.DataFrame(
        {
            "DP_raw": dp_raw,
            "OP_raw": op_raw,
            "SS_neighbors": m_cnt,
        },
        index=buildings.index,
    )
