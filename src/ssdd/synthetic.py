"""Synthetic building-footprint generators for tests and sensitivity work.

Every generator returns a :class:`geopandas.GeoDataFrame` in projected meters
(EPSG:32611 by default) so it can be passed straight into
:func:`ssdd.pipeline.compute_raw_metrics` without any reprojection.

Most geometries are simple rectangles/grids chosen so that ``SD`` and ``SS`` can
be hand-computed and the pipeline's outputs verified. Two richer generators
exercise the metric knobs that uniform rectangles cannot: ``sized_cluster``
(heterogeneous size / orientation, so the SD ``weight`` choices diverge) and
``lshape`` (a concave footprint, so ``phi_deg`` and wall-to-wall ``SS`` see a
non-rectangular outline). All also serve as parametric scaffolds for sensitivity
studies: vary ``spacing``, ``orientation_offset_deg``, ``pitch``, ``area_range``,
etc. and watch how the metrics respond.

A rectangle of size ``(w, h)`` is built with its long axis along **y** when
``h > w`` (the default). Rotation, when requested, is around the rectangle's
centroid.
"""

from __future__ import annotations

import math
from typing import Iterable

import geopandas as gpd
import numpy as np
from shapely.affinity import rotate, translate
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

    Useful baseline: both neighbor-based metrics should be 0 (SD has no
    neighbors; SS treats it as isolated).
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
    buildings share the same orientation (phi_deg identical across the grid).
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

    Saturates the separation metric: ``SS = 1 / epsilon``.
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
    area_range: tuple[float, float] | None = None,
    angle_range: tuple[float, float] | None = None,
    seed: int = 0,
    crs: int = DEFAULT_CRS,
) -> gpd.GeoDataFrame:
    """``n`` non-overlapping rectangles placed uniformly in a square.

    Useful for benchmarking and sensitivity sweeps. Buildings may be close but
    are guaranteed not to overlap via simple rejection sampling.

    ``area_range`` — if given ``(a_min, a_max)``, each footprint's area is drawn
    uniformly from that range and its ``width``/``length`` scaled to match,
    preserving the base aspect ratio. This is what makes the SD weighting choices
    (``unit`` / ``area`` / ``root_area``) diverge — on uniform-size inputs they
    differ only by a global constant.

    ``angle_range`` — if given ``(lo, hi)`` degrees, each footprint is rotated by
    a uniform draw, so orientation differences (and hence orientation-weighted SS)
    are non-degenerate.
    """
    rng = np.random.default_rng(seed)
    base_area = width * length
    geoms: list[Polygon] = []
    max_pad = (max(width, length) / 2.0 + 0.1)
    if area_range is not None:
        max_pad *= math.sqrt(area_range[1] / base_area)
    attempts = 0
    max_attempts = 50 * n
    while len(geoms) < n and attempts < max_attempts:
        attempts += 1
        w, ln = width, length
        if area_range is not None:
            s = math.sqrt(rng.uniform(*area_range) / base_area)
            w, ln = width * s, length * s
        ang = rng.uniform(*angle_range) if angle_range is not None else 0.0
        cx = rng.uniform(-extent / 2 + max_pad, extent / 2 - max_pad)
        cy = rng.uniform(-extent / 2 + max_pad, extent / 2 - max_pad)
        candidate = _rectangle(cx, cy, w, ln, angle_deg=ang)
        if any(candidate.intersects(g) for g in geoms):
            continue
        geoms.append(candidate)
    if len(geoms) < n:
        raise RuntimeError(
            f"random_cloud: only placed {len(geoms)}/{n} rectangles after "
            f"{attempts} attempts — extent {extent} is too small for the density."
        )
    return _as_gdf(geoms, crs=crs)


def sized_cluster(
    specs: Iterable[tuple[float, ...]],
    crs: int = DEFAULT_CRS,
) -> gpd.GeoDataFrame:
    """Rectangles from explicit specs — a deterministic, hand-checkable fixture
    for heterogeneous size / orientation.

    Each spec is ``(cx, cy, width, length)`` or ``(cx, cy, width, length, angle_deg)``.
    Building 0 is a convenient "focal" whose SD/SS can be computed by hand from
    the neighbour sizes and distances.
    """
    geoms = []
    for s in specs:
        cx, cy, w, ln = s[0], s[1], s[2], s[3]
        ang = s[4] if len(s) > 4 else 0.0
        geoms.append(_rectangle(cx, cy, w, ln, angle_deg=ang))
    return _as_gdf(geoms, crs=crs)


def lshape(
    arm_length: float = 30.0,
    arm_width: float = 12.0,
    angle_deg: float = 0.0,
    origin: tuple[float, float] = (0.0, 0.0),
    crs: int = DEFAULT_CRS,
) -> gpd.GeoDataFrame:
    """A single L-shaped (concave) footprint, centred on its centroid at ``origin``.

    Two arms of width ``arm_width`` and length ``arm_length`` meet at a right
    angle. Area = ``arm_width * (2 * arm_length - arm_width)``. Rotating by
    ``angle_deg`` rotates the whole footprint about ``origin``.

    Exercises the non-rectangular path: ``phi_deg`` comes from the minimum
    rotated rectangle (not a bounding box aligned to the arms), and SS measures
    true polygon-to-polygon distance to the concave outline.
    """
    poly = Polygon([
        (0.0, 0.0), (arm_length, 0.0), (arm_length, arm_width),
        (arm_width, arm_width), (arm_width, arm_length), (0.0, arm_length),
    ])
    c = poly.centroid
    poly = translate(poly, origin[0] - c.x, origin[1] - c.y)
    if angle_deg != 0.0:
        poly = rotate(poly, angle_deg, origin=(origin[0], origin[1]))
    return _as_gdf([poly], crs=crs)
