# %% [markdown]
# # 01 · Primitives
#
# Build the per-focal neighbor lists that every downstream metric reads.
# For each DINS-labeled building in Eaton + Palisades, find all neighbors within
# a centroid search radius and record the per-pair geometry:
#
#   d_wall    wall-to-wall polygon distance   (SS separation)
#   d_cent    centroid-to-centroid distance   (SD density)
#   phi_diff  folded orientation diff [0,90]  (SS orientation weight)
#   n_area    neighbor footprint area         (SD area / root-area weighting)
#
# Emits `focal.parquet` (one row per focal) and `pairs.parquet` (one row per
# focal→neighbor pair). Runs top-to-bottom as a script or cell-by-cell in Positron.

# %%
from __future__ import annotations

import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
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
except NameError:                       # Positron interactive: no __file__
    ROOT = _find_root(Path.cwd())

OUT_DIR = ROOT / "dev/py/metric-specification/_out"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FIRES = ("eaton", "palisades", "mountain")   # glass dropped: 17% missing + backwards density (gpkg kept on disk)
SEARCH_R = 250.0        # centroid meters — covers r_D up to 200 with margin

# Incumbent production metrics carried straight off the processed gpkgs so the
# 4→2 collapse can be checked against the real 4-metric set (not a proxy).
INCUMBENT = ["KD_raw", "BA_raw", "DP_raw", "OP_raw"]

print(f"root       = {ROOT}")
print(f"out_dir    = {OUT_DIR}")


def angle_difference_deg(a: float, b: float) -> float:
    """Absolute orientation difference, folded to [0, 90] deg.

    Inlined (mirrors ssdd.geometry.angle_difference_deg) so this pipeline has no
    code dependency on the package — it consumes only the processed gpkgs and
    survives package refactors.
    """
    diff = abs(a - b) % 180.0
    diff = min(diff, 180.0 - diff)
    return min(diff, 90.0)


# %%
def build_neighbor_lists(fire: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    gpkg = ROOT / f"_data/processed/{fire}/{fire}_buildings.gpkg"
    bld = gpd.read_file(gpkg)
    if "DAMAGE" not in bld.columns:
        raise RuntimeError(f"{gpkg} has no DAMAGE column — join with DINS first")

    bld = bld[bld["DAMAGE"].notna()].copy().reset_index(drop=True)
    bld["destroyed"] = (bld["DAMAGE"] == "Destroyed (>50%)").astype(int)
    bld["fire"] = fire
    if bld.crs is None or bld.crs.to_epsg() != 32611:
        bld = bld.to_crs(epsg=32611)
    print(f"[{fire}] {len(bld):,} DINS-labeled buildings "
          f"({bld['destroyed'].sum():,} destroyed)")

    polys = bld.geometry.values
    cx = bld["cent_x"].to_numpy(float)
    cy = bld["cent_y"].to_numpy(float)
    phi = bld["phi_deg"].to_numpy(float)
    area = bld["bld_area"].to_numpy(float)
    ids = bld["ssdd_id"].to_numpy()

    cent_pts = [Point(x, y) for x, y in zip(cx, cy)]
    tree = STRtree(cent_pts)

    f_id, n_id, dw, dc, dphi, na = [], [], [], [], [], []
    for i in range(len(bld)):
        for j in tree.query(cent_pts[i].buffer(SEARCH_R)):
            if j == i:
                continue
            d_c = float(np.hypot(cx[i] - cx[j], cy[i] - cy[j]))
            if d_c > SEARCH_R:
                continue
            f_id.append(ids[i]); n_id.append(ids[j])
            dw.append(float(polys[i].distance(polys[j]))); dc.append(d_c)
            dphi.append(angle_difference_deg(float(phi[i]), float(phi[j])))
            na.append(area[j])
        if i and i % 3000 == 0:
            print(f"  {i:,}/{len(bld):,} focals — {len(f_id):,} pairs")

    pairs = pd.DataFrame({"focal_id": f_id, "neighbor_id": n_id, "fire": fire,
                          "d_wall": dw, "d_cent": dc, "phi_diff": dphi, "n_area": na})
    keep = ["ssdd_id", "fire", "destroyed", "bld_area", "phi_deg", "cent_x", "cent_y"]
    keep += [c for c in INCUMBENT if c in bld.columns]
    focal = bld[keep].rename(columns={"ssdd_id": "focal_id"})
    return focal, pairs


# %%
focals, all_pairs = [], []
for fire in FIRES:
    f, p = build_neighbor_lists(fire)
    focals.append(f); all_pairs.append(p)

focal_df = pd.concat(focals, ignore_index=True)
pairs_df = pd.concat(all_pairs, ignore_index=True)

focal_df.to_parquet(OUT_DIR / "focal.parquet", index=False)
pairs_df.to_parquet(OUT_DIR / "pairs.parquet", index=False)
print(f"\nWrote focal.parquet ({len(focal_df):,} rows) and "
      f"pairs.parquet ({len(pairs_df):,} rows)")
