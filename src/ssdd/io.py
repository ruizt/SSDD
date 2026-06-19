"""I/O helpers: reading buildings/DINS, CRS handling, spatial join."""

from __future__ import annotations

from typing import Literal, Optional, cast

import geopandas as gpd
import numpy as np

JoinHow = Literal["left", "inner", "right"]


def read_buildings(path: str, layer: Optional[str] = None) -> gpd.GeoDataFrame:
    """Read building footprints, drop null geometries, repair with buffer(0).

    Keeps only Polygon / MultiPolygon features.
    """
    gdf = gpd.read_file(path, layer=layer) if layer else gpd.read_file(path)
    gdf = cast(gpd.GeoDataFrame, gdf[gdf.geometry.notnull()].copy())
    gdf = cast(gpd.GeoDataFrame,
               gdf[gdf.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy())
    gdf["geometry"] = gdf["geometry"].buffer(0)
    return cast(gpd.GeoDataFrame, gdf.reset_index(drop=True))


def read_dins(path: str, layer: Optional[str] = None) -> gpd.GeoDataFrame:
    """Read DINS structure points. Keeps only Point geometries."""
    gdf = gpd.read_file(path, layer=layer) if layer else gpd.read_file(path)
    gdf = cast(gpd.GeoDataFrame, gdf[gdf.geometry.notnull()].copy())
    gdf = cast(gpd.GeoDataFrame,
               gdf[gdf.geometry.geom_type == "Point"].copy())
    return cast(gpd.GeoDataFrame, gdf.reset_index(drop=True))


def ensure_projected_meters(gdf: gpd.GeoDataFrame, target_epsg: int) -> gpd.GeoDataFrame:
    """Reproject to ``target_epsg`` if needed. Raises if input has no CRS."""
    if gdf.crs is None:
        raise ValueError(
            "Input has no CRS defined. Set it in QGIS (Layer Properties → Source) first."
        )
    if gdf.crs.to_epsg() != target_epsg:
        gdf = gdf.to_crs(epsg=target_epsg)
    return gdf


def join_dins(
    buildings: gpd.GeoDataFrame,
    dins: gpd.GeoDataFrame,
    keep_cols: Optional[list[str]] = None,
    how: JoinHow = "left",
    predicate: str = "intersects",
    id_col: str = "ssdd_id",
) -> gpd.GeoDataFrame:
    """Attach DINS attributes to buildings via point-in-polygon spatial join.

    Parameters
    ----------
    buildings
        Building footprints. Must already carry an identifier column ``id_col``
        (use :func:`ssdd.pipeline.compute_raw_metrics` upstream, or set one
        explicitly).
    dins
        DINS points. CRS must match ``buildings`` — caller should reproject first.
    keep_cols
        Optional whitelist of DINS attribute columns to keep. If ``None``, all
        non-geometry DINS columns are carried over.
    how
        ``"left"`` keeps all buildings, NaN where no DINS point falls inside.
        ``"inner"`` keeps only DINS-overlapping buildings (the "burned" subset).
    predicate
        Spatial predicate. Default ``"intersects"`` matches the R workflow's
        ``st_join`` semantics for point-in-polygon.
    id_col
        Column on ``buildings`` used to de-duplicate when multiple DINS points
        fall in the same footprint (rare; keeps the first).

    Returns
    -------
    GeoDataFrame
        ``buildings`` with DINS attributes appended. The ``index_right`` column
        from sjoin is dropped.
    """
    if buildings.crs != dins.crs:
        raise ValueError(
            f"CRS mismatch: buildings={buildings.crs}, dins={dins.crs}. "
            "Reproject one to match before joining."
        )
    if id_col not in buildings.columns:
        raise KeyError(
            f"Buildings is missing identifier column {id_col!r}. "
            "Call compute_raw_metrics first or assign one yourself."
        )

    if keep_cols is not None:
        dins_subset = dins[keep_cols + ["geometry"]].copy()
    else:
        dins_subset = dins.copy()

    joined = gpd.sjoin(buildings, cast(gpd.GeoDataFrame, dins_subset),
                       how=how, predicate=predicate)
    joined = joined.drop(columns=["index_right"], errors="ignore")
    joined = joined.drop_duplicates(subset=[id_col], keep="first").reset_index(drop=True)
    return cast(gpd.GeoDataFrame, joined)


def add_building_id(buildings: gpd.GeoDataFrame, id_col: str = "ssdd_id") -> gpd.GeoDataFrame:
    """Ensure ``buildings`` has a stable integer identifier column.

    Default name ``ssdd_id`` avoids case-insensitive collisions with common
    source-data identifiers like LARIAC6's ``BLD_ID`` when writing to
    GeoPackage / SQLite, which fold column names case-insensitively.
    """
    if id_col not in buildings.columns:
        buildings = buildings.copy()
        buildings[id_col] = np.arange(len(buildings), dtype=np.int64)
    return buildings
