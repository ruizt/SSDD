# %% [markdown]
# # 00 · Prepare processed gpkgs  (package-free)
#
# Turns a raw footprint file + DINS points into the processed
# `_data/processed/{fire}/{fire}_buildings.gpkg` that `01_primitives` consumes —
# WITHOUT importing the `ssdd` package. Computes exactly the fields the metric
# grids need:
#
#   ssdd_id · bld_area · phi_deg · cent_x · cent_y · DAMAGE (point-in-polygon join)
#
# `phi_deg` reproduces `ssdd.geometry.dominant_orientation_degrees` (longest edge
# of the minimum rotated rectangle) so it matches the Eaton/Palisades gpkgs bit
# for bit. Extend FIRES to prep new fires. Eaton/Palisades keep their existing
# (package-made) gpkgs; pass their LARIAC shapefiles here if you want to
# regenerate them the same package-free way.

# %%
from __future__ import annotations

import math
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.strtree import STRtree

warnings.filterwarnings("ignore")


def _find_root(start: Path) -> Path:
    p = start.resolve()
    while p != p.parent:
        if (p / "_data").exists() and (p / "dev").exists():
            return p
        p = p.parent
    return start.resolve()


try:
    ROOT = _find_root(Path(__file__).parent)
except NameError:
    ROOT = _find_root(Path.cwd())

EPSG = 32611

# fire -> (raw footprints, DINS points)
FIRES = {
    "glass":    ("_data/raw/buildings/MSFootprints_2020_Glass.geojson",
                 "_data/raw/dins/DINS_2020_Glass.geojson"),
    "mountain": ("_data/raw/buildings/USAStructures_2024_Mountain.geojson",
                 "_data/raw/dins/DINS_2024_Mountain.geojson"),
}


def dominant_orientation_degrees(poly) -> float:
    """Longest-edge azimuth of the minimum rotated rectangle, in [0,180).
    Copy of ssdd.geometry.dominant_orientation_degrees (kept inline for
    package independence)."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, module=r"shapely\..*")
        mrr = poly.minimum_rotated_rectangle
    c = list(mrr.exterior.coords)
    best = max(((math.hypot(c[k+1][0]-c[k][0], c[k+1][1]-c[k][1]),
                 c[k+1][0]-c[k][0], c[k+1][1]-c[k][1]) for k in range(4)),
               key=lambda t: t[0])
    return math.degrees(math.atan2(best[2], best[1])) % 180.0


def prepare(fire: str, bld_path: str, dins_path: str) -> None:
    bld = gpd.read_file(ROOT / bld_path)
    if bld.crs is None or bld.crs.to_epsg() != EPSG:
        bld = bld.to_crs(EPSG)
    bld = bld[bld.geometry.notna() & ~bld.geometry.is_empty].reset_index(drop=True)

    bld["ssdd_id"] = np.arange(len(bld))
    bld["bld_area"] = bld.geometry.area
    cent = bld.geometry.centroid
    bld["cent_x"] = cent.x.to_numpy()
    bld["cent_y"] = cent.y.to_numpy()
    bld["phi_deg"] = [dominant_orientation_degrees(g) for g in bld.geometry.values]

    dins = gpd.read_file(ROOT / dins_path)
    if dins.crs is None or dins.crs.to_epsg() != EPSG:
        dins = dins.to_crs(EPSG)
    # point-in-polygon join (matches ssdd.io.join_dins: predicate="intersects", left)
    dmg = np.full(len(bld), None, dtype=object)
    polys = list(bld.geometry.values)
    tree = STRtree(polys)
    for pt, dm in zip(dins.geometry.values, dins["DAMAGE"].values):
        for i in tree.query(pt):
            if polys[i].contains(pt) or polys[i].intersects(pt):
                if dmg[i] is None:
                    dmg[i] = dm
                break
    bld["DAMAGE"] = dmg

    out = ROOT / f"_data/processed/{fire}/{fire}_buildings.gpkg"
    out.parent.mkdir(parents=True, exist_ok=True)
    keep = ["ssdd_id", "bld_area", "phi_deg", "cent_x", "cent_y", "DAMAGE", "geometry"]
    bld[keep].to_file(out, layer="buildings_raw", driver="GPKG")
    n_lab = int(bld["DAMAGE"].notna().sum())
    print(f"[{fire}] {len(bld):,} buildings, {n_lab:,} DINS-labeled -> {out.name}")


# %%
for fire, (bp, dp) in FIRES.items():
    prepare(fire, bp, dp)
