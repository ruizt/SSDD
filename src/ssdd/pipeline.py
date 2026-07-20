"""End-to-end orchestration: buildings GeoDataFrame -> raw SSDD metrics.

``compute_raw_metrics`` attaches per-building attributes (id, area, dominant
orientation, centroid coords) and the two raw metrics ``SD`` and
``SS``. Normalization, blending and modeling are downstream choices.
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
    compute_SD_series,
    compute_SS_series,
)


@dataclass
class RawMetricParams:
    """Parameters of the two raw metrics.

    Defaults follow the metric-specification experiments in
    ``dev/py/metric-specification/`` (Eaton, Palisades, Mountain;
    structure-level AUC against DINS damage):

    * ``r_D = 200`` — separation improves monotonically across the swept grid
      (pooled AUC 0.646 @50 → 0.706 @200) and is the per-fire optimum for Eaton
      and Mountain (Palisades prefers 100). NOTE: 200 m is the top of the swept
      range, so larger radii remain untested.
    * ``kernel = "uniform"`` — kernel shape is immaterial (all forms within
      ~0.02 AUC; the per-fire winner varies), so the simplest is preferred.
    * ``weight = "root_area"`` — best pooled separation (0.706 vs 0.697 unit,
      0.670 area at r_D = 200) and never worse than 2nd in any single fire,
      while ``unit`` is worst in Eaton and ``area`` is worst in Palisades.
    * ``r_S = 50`` — the nearest-neighbour separation form is r_S-invariant
      (AUC 0.619 at r_S = 25 / 50 / 75); 50 is retained as a sensible middle.
    * ``epsilon = 0.5`` — **not** established by these experiments: held fixed
      throughout and never swept. Bounds SS at 1/epsilon for touching walls.
    * ``agg = "nn"`` / ``orient = "flat"`` — the plain nearest-neighbour form
      dominates every averaged / orientation-weighted alternative (pooled AUC
      0.619 vs <= 0.59). The alternative aggregations and orientation weights
      remain available as parameters (the legacy DP is ``agg="power1"``, the
      legacy OP is ``agg="power1", orient="gauss"``); ``sigma`` only
      applies when ``orient="gauss"``.
    """

    r_D: float = 200.0
    r_S: float = 50.0
    epsilon: float = 0.5
    kernel: str = "uniform"
    weight: str = "root_area"
    agg: str = "nn"
    orient: str = "flat"
    sigma: float = 10.0


def compute_raw_metrics(
    buildings: gpd.GeoDataFrame,
    params: RawMetricParams | None = None,
    id_col: str = "ssdd_id",
    progress: bool = True,
) -> gpd.GeoDataFrame:
    """Compute per-building attributes and the raw SD / SS metrics.

    Parameters
    ----------
    buildings
        Footprint polygons in a projected meters CRS
        (:func:`ssdd.io.ensure_projected_meters` upstream).
    params
        :class:`RawMetricParams`; defaults follow the metric-specification
        experiments.
    id_col
        Name of the stable integer id column added if missing.
    progress
        Show tqdm progress bars.

    Returns
    -------
    GeoDataFrame
        Input plus ``ssdd_id`` (if missing), ``bld_area``, ``phi_deg``
        (dominant footprint orientation — an attribute for downstream analyses,
        not a metric input), ``cent_x`` / ``cent_y`` (true polygon centroid
        coords in the input CRS), ``SD`` and ``SS``.
    """
    p = params or RawMetricParams()

    bld = add_building_id(buildings, id_col=id_col).copy()
    bld["bld_area"] = bld.geometry.area

    if progress:
        tqdm.pandas(desc="  orientation")
        bld["phi_deg"] = bld.geometry.progress_apply(dominant_orientation_degrees)
    else:
        bld["phi_deg"] = bld.geometry.apply(dominant_orientation_degrees)

    polys = bld.geometry.values
    cents = bld.geometry.centroid.values
    tree_polys = STRtree(polys)
    tree_cents = STRtree(cents)

    bld["cent_x"] = [c.x for c in cents]
    bld["cent_y"] = [c.y for c in cents]

    bld["SD"] = compute_SD_series(
        bld,
        r_D=p.r_D,
        kernel=p.kernel,
        weight=p.weight,
        tree_cents=tree_cents,
        cents=cents,
        areas=bld["bld_area"].to_numpy() if p.weight != "unit" else None,
        progress=progress,
    )
    bld["SS"] = compute_SS_series(
        bld,
        r_S=p.r_S,
        epsilon=p.epsilon,
        agg=p.agg,
        orient=p.orient,
        sigma=p.sigma,
        phi_deg=bld["phi_deg"].to_numpy(),
        tree_polys=tree_polys,
        polys=polys,
        progress=progress,
    )

    return cast(gpd.GeoDataFrame, bld)
