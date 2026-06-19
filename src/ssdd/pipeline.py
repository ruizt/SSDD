"""End-to-end raw-metric computation.

This is the function you'd typically call from a notebook or wrapper script:
give it a building GeoDataFrame (already in projected meters) and parameters,
get back a copy with all four raw metrics + supporting fields.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import geopandas as gpd
import pandas as pd
from shapely.strtree import STRtree
from tqdm.auto import tqdm

from .geometry import dominant_orientation_degrees
from .io import add_building_id
from .metrics import (
    compute_BA_series,
    compute_KD_series,
    compute_NN_proximity,
    compute_SS_terms_df,
)


@dataclass
class RawMetricParams:
    """Spatial parameters that define neighborhood size and shape."""

    r_D: float = 100.0
    r_S: float = 50.0
    epsilon: float = 0.5
    sigma_theta: float = 15.0
    kernel: str = "quartic"
    weight_by_area: bool = False
    r_NN: float = 200.0


def compute_raw_metrics(
    buildings: gpd.GeoDataFrame,
    params: RawMetricParams | None = None,
    id_col: str = "ssdd_id",
    progress: bool = True,
    **kwargs,
) -> gpd.GeoDataFrame:
    """Compute KD_raw, BA_raw, DP_raw, OP_raw and supporting fields.

    Parameters
    ----------
    buildings
        GeoDataFrame of building polygons in a projected CRS with meter units.
    params
        :class:`RawMetricParams`. Individual fields can also be passed as
        keyword arguments and will override the dataclass values.
    id_col
        Identifier column to ensure on the output. Created from ``arange`` if
        missing.

    Returns
    -------
    GeoDataFrame
        Copy of ``buildings`` with these columns added:

        ``ssdd_id`` (if missing), ``bld_area``, ``phi_deg``,
        ``cent_x``, ``cent_y`` (representative-point coords in the input CRS),
        ``KD_raw``, ``BA_raw``, ``DP_raw``, ``OP_raw``, ``SS_neighbors``,
        ``dist_to_nearest_building``, ``bearing_to_nearest_building``.
    """
    p = params or RawMetricParams()
    if kwargs:
        # Allow per-call overrides without forcing the caller to build a dataclass.
        p = RawMetricParams(**{**p.__dict__, **kwargs})

    bld = add_building_id(buildings, id_col=id_col).copy()
    bld["bld_area"] = bld.geometry.area

    if progress:
        tqdm.pandas(desc="  orientation")
        bld["phi_deg"] = bld.geometry.progress_apply(dominant_orientation_degrees)
    else:
        bld["phi_deg"] = bld.geometry.apply(dominant_orientation_degrees)

    polys = bld.geometry.values
    rep_pts = bld.geometry.representative_point().values
    tree_polys = STRtree(polys)
    tree_pts = STRtree(rep_pts)

    bld["cent_x"] = [p.x for p in rep_pts]
    bld["cent_y"] = [p.y for p in rep_pts]

    bld["KD_raw"] = compute_KD_series(
        bld,
        r_D=p.r_D,
        kernel=p.kernel,
        weight_by_area=p.weight_by_area,
        tree_pts=tree_pts,
        rep_pts=rep_pts,
        areas=bld["bld_area"].to_numpy() if p.weight_by_area else None,
        progress=progress,
    )
    bld["BA_raw"] = compute_BA_series(
        bld,
        r_D=p.r_D,
        tree_polys=tree_polys,
        polys=polys,
        progress=progress,
    )
    ss = compute_SS_terms_df(
        bld,
        r_S=p.r_S,
        epsilon=p.epsilon,
        sigma_theta=p.sigma_theta,
        phi_deg=cast(pd.Series, bld["phi_deg"]),
        tree_polys=tree_polys,
        polys=polys,
        progress=progress,
    )
    bld[["DP_raw", "OP_raw", "SS_neighbors"]] = ss

    nn = compute_NN_proximity(
        bld,
        r_NN=p.r_NN,
        tree_polys=tree_polys,
        polys=polys,
        rep_pts=rep_pts,
        progress=progress,
    )
    bld[["dist_to_nearest_building", "bearing_to_nearest_building"]] = nn

    return bld
