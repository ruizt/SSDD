"""Synthetic building-footprint generators for tests and sensitivity work.

Every generator returns a :class:`geopandas.GeoDataFrame` in projected meters
(EPSG:32611 by default) so it can be passed straight into
:func:`ssdd.pipeline.compute_raw_metrics` without any reprojection.

The geometries are simple — rectangles, grids, clusters — chosen so that
``KD_raw``, ``BA_raw``, ``DP_raw`` and ``OP_raw`` can be hand-computed and the
pipeline's outputs can be verified. They also serve as parametric scaffolds for
sensitivity studies: vary ``spacing``, ``orientation_offset_deg``, ``pitch``,
etc. and watch how raw and blended metrics respond.

A rectangle of size ``(w, h)`` is built with its long axis along **y** when
``h > w`` (the default). Rotation, when requested, is around the rectangle's
centroid.
"""

from __future__ import annotations

from typing import Iterable

import geopandas as gpd
import numpy as np
from shapely.affinity import rotate
from shapely.geometry import Polygon

DEFAULT_CRS = 32611  # UTM 11N, meters — matches the SSDD analysis CRS.


def _rectangle(cx: float, cy: float, w: float, h: float, angle_deg: float = 0.0) -> Polygon:
    """Axis-aligned rectangle centered at (cx, cy), optionally rotated about its centroid."""
    poly = Polygon([
        (cx - w / 2, cy - h / 2),
        (cx + w / 2, cy - h / 2),
        (cx + w / 2, cy + h / 2),
        (cx - w / 2, cy + h / 2),
    ])
    if angle_deg != 0.0:
        poly = rotate(poly, angle_deg, origin=(cx, cy))
    return poly


def _as_gdf(geoms: Iterable[Polygon], crs: int = DEFAULT_CRS) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame({"geometry": list(geoms)}, geometry="geometry", crs=f"EPSG:{crs}")


def isolated_building(
    width: float = 10.0,
    length: float = 20.0,
    angle_deg: float = 0.0,
    origin: tuple[float, float] = (0.0, 0.0),
    crs: int = DEFAULT_CRS,
) -> gpd.GeoDataFrame:
    """A single rectangle, no neighbors.

    Useful baseline: every neighbor-based metric should be 0; ``BA_raw`` should
    equal the building's own area divided by the Minkowski-sum window area.
    """
    return _as_gdf(
        [_rectangle(origin[0], origin[1], width, length, angle_deg)],
        crs=crs,
    )


def pair(
    spacing: float = 10.0,
    orientation_offset_deg: float = 0.0,
    width: float = 10.0,
    length: float = 20.0,
    crs: int = DEFAULT_CRS,
) -> gpd.GeoDataFrame:
    """Two rectangles separated by ``spacing`` (pre-rotation wall-to-wall, meters).

    B1 sits left of origin, B2 sits right of origin. B2 is rotated by
    ``orientation_offset_deg`` about its own centroid. When the offset is 0,
    the wall-to-wall distance equals ``spacing`` exactly; with non-zero
    rotation, the wall-to-wall distance is smaller because the rotated
    rectangle's corner extends inward. Use :func:`shapely.geometry.distance`
    on the returned geometries to get the actual wall-to-wall distance for
    expected-value calculations in tests.
    """
    cx1 = -(width + spacing) / 2.0
    cx2 = +(width + spacing) / 2.0
    b1 = _rectangle(cx1, 0.0, width, length, angle_deg=0.0)
    b2 = _rectangle(cx2, 0.0, width, length, angle_deg=orientation_offset_deg)
    return _as_gdf([b1, b2], crs=crs)


def grid(
    n: int = 3,
    pitch: float = 20.0,
    width: float = 10.0,
    length: float = 20.0,
    angle_deg: float = 0.0,
    crs: int = DEFAULT_CRS,
) -> gpd.GeoDataFrame:
    """An ``n x n`` regular grid of rectangles, centers on a square lattice.

    ``pitch`` is the center-to-center distance (not wall-to-wall). All
    buildings share the same orientation, so orientation differences are 0
    and OP equals DP.
    """
    geoms = []
    offset = (n - 1) * pitch / 2.0
    for i in range(n):
        for j in range(n):
            cx = i * pitch - offset
            cy = j * pitch - offset
            geoms.append(_rectangle(cx, cy, width, length, angle_deg))
    return _as_gdf(geoms, crs=crs)


def touching_pair(
    width: float = 10.0,
    length: float = 20.0,
    orientation_offset_deg: float = 0.0,
    crs: int = DEFAULT_CRS,
) -> gpd.GeoDataFrame:
    """Two rectangles sharing a wall (wall-to-wall = 0 m).

    Saturates the inverse-distance proxies: ``DP_raw = 1 / epsilon`` and
    ``OP_raw = g(theta) / epsilon``.
    """
    return pair(
        spacing=0.0,
        orientation_offset_deg=orientation_offset_deg,
        width=width,
        length=length,
        crs=crs,
    )


def random_cloud(
    n: int = 100,
    extent: float = 500.0,
    width: float = 10.0,
    length: float = 20.0,
    seed: int = 0,
    crs: int = DEFAULT_CRS,
) -> gpd.GeoDataFrame:
    """``n`` non-overlapping axis-aligned rectangles placed uniformly in a square.

    Useful for benchmarking and sensitivity sweeps. Buildings may be close but
    are guaranteed not to overlap via simple rejection sampling.
    """
    rng = np.random.default_rng(seed)
    geoms: list[Polygon] = []
    pad_x = width / 2.0 + 0.1
    pad_y = length / 2.0 + 0.1
    attempts = 0
    max_attempts = 50 * n
    while len(geoms) < n and attempts < max_attempts:
        attempts += 1
        cx = rng.uniform(-extent / 2 + pad_x, extent / 2 - pad_x)
        cy = rng.uniform(-extent / 2 + pad_y, extent / 2 - pad_y)
        candidate = _rectangle(cx, cy, width, length)
        if any(candidate.intersects(g) for g in geoms):
            continue
        geoms.append(candidate)
    if len(geoms) < n:
        raise RuntimeError(
            f"random_cloud: only placed {len(geoms)}/{n} rectangles after "
            f"{attempts} attempts — extent {extent} is too small for the density."
        )
    return _as_gdf(geoms, crs=crs)
