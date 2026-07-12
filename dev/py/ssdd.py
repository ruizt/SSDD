#!/usr/bin/env python3
"""
Structure Separation Distance Density (SSDD) — historical prototype
====================================================================
This single-file script is the **original prototype** that defined the SSDD
metric, including its convex-combination blending into SD, SS, and the final
SSDD score. It was derived from the exploratory notebook at
``dev/notebooks/SSDD.ipynb`` (a notebook → script export).

It is retained for **reference only**. The production package in ``src/ssdd/``
deliberately stops at the four raw metrics (KD, BA, DP, OP) and leaves
normalization, blending and predictive modeling to the user. If you want to
see how the originally proposed blending worked, read on; if you want to
compute metrics for new data, use ``src/ssdd_compute.py`` or the package API.

Two components:
  SD (Structure Density)  -- kernel density + basal area fraction, blended via alpha_D
  SS (Separation)         -- inverse-distance + orientation-weighted inverse-distance, blended via alpha_S
  SSDD                    -- linear blend of SD and SS via beta

Usage (CLI):
    python ssdd.py --input /path/to/Buildings.shp --output /path/to/output --run-name my_run

Usage (edit parameters below):
    Set INPUT_PATH, OUT_DIR, RUN_NAME etc. in the PARAMETERS section, then:
    python ssdd.py

Dependencies:
    geopandas, shapely, numpy, pandas, matplotlib, tqdm
"""

import os
import math
import time
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
from shapely.strtree import STRtree


# =========================
# PARAMETERS
# Edit these directly, or override via CLI arguments (see --help).
# =========================

INPUT_PATH   = r"/path/to/Buildings.shp"   # SHP / GPKG / GeoJSON / FGDB
INPUT_LAYER  = None                          # e.g. "buildings" for GPKG; None for SHP

OUT_DIR      = r"/path/to/output"
RUN_NAME     = "my_run"

TARGET_EPSG  = 3310   # NAD83 / California Albers (equal-area, meters)

# ---------- SD knobs ----------
r_D            = 100.0      # SD buffer radius (m)
alpha_D        = 0.5        # 1 = kernel density only, 0 = basal area only
kernel_type    = "quartic"
weight_by_area = False      # True => weight neighbor contributions by footprint area

# ---------- SS knobs ----------
r_S         = 50.0    # SS neighbor search radius (m)
epsilon     = 0.5     # distance floor to avoid division-by-zero (m)
sigma_theta = 15.0    # orientation tolerance (deg); smaller => orientation matters more
alpha_S     = 0.5     # 1 = distance only, 0 = orientation-weighted distance only

# ---------- Master blend ----------
beta = 0.5   # 1 = SD only, 0 = SS only

# ---------- Normalization ----------
NORM_METHOD = "robust"   # "minmax" or "robust"
P_LOW, P_HIGH = 2, 98

# ---------- QA ----------
MAKE_QA_PLOTS = True


# =========================
# CLI ARGUMENT PARSING
# =========================

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute SSDD metrics for building footprints.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",       default=None, help="Input shapefile/GPKG/GeoJSON path")
    p.add_argument("--layer",       default=None, help="Layer name (GPKG/FGDB); omit for SHP")
    p.add_argument("--output",      default=None, help="Output directory")
    p.add_argument("--run-name",    default=None, help="Run name prefix for output files")
    p.add_argument("--epsg",        type=int,   default=None, help="Target EPSG code")
    p.add_argument("--r-d",         type=float, default=None, help="SD buffer radius (m)")
    p.add_argument("--alpha-d",     type=float, default=None, help="SD blend alpha (0–1)")
    p.add_argument("--r-s",         type=float, default=None, help="SS search radius (m)")
    p.add_argument("--epsilon",     type=float, default=None, help="Distance floor (m)")
    p.add_argument("--sigma-theta", type=float, default=None, help="Orientation tolerance (deg)")
    p.add_argument("--alpha-s",     type=float, default=None, help="SS blend alpha (0–1)")
    p.add_argument("--beta",        type=float, default=None, help="SD/SS master blend (0–1)")
    p.add_argument("--norm",        default=None, choices=["minmax", "robust"],
                   help="Normalization method")
    p.add_argument("--no-qa-plots", action="store_true", help="Suppress QA histogram plots")
    return p.parse_args()


# =========================
# UTILITY FUNCTIONS
# =========================

def read_buildings(path: str, layer):
    if layer is None:
        return gpd.read_file(path)
    return gpd.read_file(path, layer=layer)


def ensure_projected_meters(gdf: gpd.GeoDataFrame, target_epsg: int) -> gpd.GeoDataFrame:
    if gdf.crs is None:
        raise ValueError(
            "Input has no CRS defined. Set it in QGIS (Layer Properties → Source) first."
        )
    if gdf.crs.to_epsg() != target_epsg:
        gdf = gdf.to_crs(epsg=target_epsg)
    return gdf


def norm_series(x: pd.Series, method="minmax", p_low=2, p_high=98) -> pd.Series:
    x = x.astype(float)
    if method == "minmax":
        lo, hi = np.nanmin(x), np.nanmax(x)
        if hi == lo:
            return pd.Series(np.zeros(len(x)), index=x.index)
        return (x - lo) / (hi - lo)
    if method == "robust":
        lo, hi = np.nanpercentile(x, p_low), np.nanpercentile(x, p_high)
        if hi == lo:
            return pd.Series(np.zeros(len(x)), index=x.index)
        x_clip = np.clip(x, lo, hi)
        return (x_clip - lo) / (hi - lo)
    raise ValueError(f"Unknown normalization method: {method}")


def quartic_kernel(u: float) -> float:
    """K(u) = (1 - u^2)^2 for u in [0, 1], else 0."""
    if u < 0 or u > 1:
        return 0.0
    return (1.0 - u * u) ** 2


def kernel_value(u: float, kernel="quartic") -> float:
    if kernel == "quartic":
        return quartic_kernel(u)
    raise ValueError(f"Unknown kernel type: {kernel}")


def dominant_orientation_degrees(poly) -> float:
    """Angle of the longest edge of the minimum rotated rectangle, in [0, 180)."""
    mrr = poly.minimum_rotated_rectangle
    coords = list(mrr.exterior.coords)
    edges = []
    for k in range(4):
        x1, y1 = coords[k]
        x2, y2 = coords[k + 1]
        dx, dy = x2 - x1, y2 - y1
        edges.append((math.hypot(dx, dy), dx, dy))
    _, dx, dy = max(edges, key=lambda t: t[0])
    return math.degrees(math.atan2(dy, dx)) % 180.0


def angle_difference_deg(a: float, b: float) -> float:
    diff = abs(a - b) % 180.0
    diff = min(diff, 180.0 - diff)
    return min(diff, 90.0)


def orientation_factor(theta_deg: float, sigma_theta: float) -> float:
    """g(theta) = exp(-(theta / sigma)^2)."""
    if sigma_theta <= 0:
        return 1.0
    return math.exp(-((theta_deg / sigma_theta) ** 2))


def tree_query_indices(tree: STRtree, query_geom):
    """Return integer indices from STRtree.query(), Shapely 2.x compatible."""
    res = tree.query(query_geom)
    if isinstance(res, (list, tuple, np.ndarray)) and len(res) > 0:
        if isinstance(res[0], (int, np.integer)):
            return np.asarray(res, dtype=int)
    idx_map = {id(g): i for i, g in enumerate(tree.geometries)}
    return np.asarray([idx_map[id(g)] for g in res], dtype=int)


def qa_hist(series: pd.Series, title: str):
    plt.figure()
    plt.hist(series.dropna().values, bins=40)
    plt.title(title)
    plt.xlabel(series.name)
    plt.ylabel("Count")
    plt.show()


# =========================
# MAIN
# =========================

def main():
    # Allow global parameters to be overridden by CLI
    global INPUT_PATH, INPUT_LAYER, OUT_DIR, RUN_NAME, TARGET_EPSG
    global r_D, alpha_D, kernel_type, weight_by_area
    global r_S, epsilon, sigma_theta, alpha_S
    global beta, NORM_METHOD, P_LOW, P_HIGH, MAKE_QA_PLOTS

    args = parse_args()
    if args.input is not None:        INPUT_PATH    = args.input
    if args.layer is not None:        INPUT_LAYER   = args.layer
    if args.output is not None:       OUT_DIR       = args.output
    if args.run_name is not None:     RUN_NAME      = args.run_name
    if args.epsg is not None:         TARGET_EPSG   = args.epsg
    if args.r_d is not None:          r_D           = args.r_d
    if args.alpha_d is not None:      alpha_D       = args.alpha_d
    if args.r_s is not None:          r_S           = args.r_s
    if args.epsilon is not None:      epsilon       = args.epsilon
    if args.sigma_theta is not None:  sigma_theta   = args.sigma_theta
    if args.alpha_s is not None:      alpha_S       = args.alpha_s
    if args.beta is not None:         beta          = args.beta
    if args.norm is not None:         NORM_METHOD   = args.norm
    if args.no_qa_plots:              MAKE_QA_PLOTS = False

    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load, clean, and reproject
    # ------------------------------------------------------------------
    print("Loading buildings...")
    bld = read_buildings(INPUT_PATH, INPUT_LAYER)
    bld = bld[bld.geometry.notnull()].copy()
    bld = bld[bld.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    bld["geometry"] = bld["geometry"].buffer(0)   # fix common geometry issues

    print(f"  Input CRS  : {bld.crs}")
    bld = ensure_projected_meters(bld, TARGET_EPSG)
    print(f"  Analysis CRS: {bld.crs}")

    if "bld_id" not in bld.columns:
        bld["bld_id"] = np.arange(len(bld))
    bld["bld_area"] = bld.geometry.area
    bld["rep_pt"]   = bld.geometry.representative_point()
    print(f"  Loaded {len(bld):,} buildings")

    # ------------------------------------------------------------------
    # 2. Dominant building orientation
    # ------------------------------------------------------------------
    print("Computing building orientations...")
    tqdm.pandas(desc="  orientation")
    bld["phi_deg"] = bld.geometry.progress_apply(dominant_orientation_degrees)

    # ------------------------------------------------------------------
    # 3. Spatial indexes
    # ------------------------------------------------------------------
    print("Building spatial indexes...")
    polys      = bld.geometry.values
    pts        = bld["rep_pt"].values
    tree_polys = STRtree(polys)
    tree_pts   = STRtree(pts)

    # ------------------------------------------------------------------
    # 4. SD: Kernel Density (KD) + Basal Area Fraction (BA)
    # ------------------------------------------------------------------

    def compute_KD(i: int) -> float:
        """
        KD_i = (1 / (pi * r_D^2)) * sum_j  w_j * K(dist(ci, cj) / r_D)
        for all neighbors j within r_D of building i's representative point.
        """
        ci   = pts[i]
        idxs = tree_query_indices(tree_pts, ci.buffer(r_D))
        total = 0.0
        for j in idxs:
            if j == i:
                continue
            dist = ci.distance(pts[j])
            if dist > r_D:
                continue
            u = dist / r_D
            w = float(bld.iloc[j]["bld_area"]) if weight_by_area else 1.0
            total += w * kernel_value(u, kernel=kernel_type)
        return total / (math.pi * r_D * r_D)

    def compute_BA(i: int) -> float:
        """
        BA_i = sum_j area(P_j ∩ buffer(P_i, r_D)) / area(buffer(P_i, r_D))
        """
        win      = polys[i].buffer(r_D)
        win_area = win.area
        if win_area <= 0:
            return 0.0
        idxs = tree_query_indices(tree_polys, win)
        inter_area_sum = 0.0
        for j in idxs:
            inter = polys[j].intersection(win)
            if not inter.is_empty:
                inter_area_sum += inter.area
        return inter_area_sum / win_area

    print("Computing SD components...")
    bld["KD_raw"] = [compute_KD(i) for i in tqdm(range(len(bld)), desc="  KD (kernel density)")]
    bld["BA_raw"] = [compute_BA(i) for i in tqdm(range(len(bld)), desc="  BA (basal area fraction)")]

    bld["KD"] = norm_series(bld["KD_raw"], method=NORM_METHOD, p_low=P_LOW, p_high=P_HIGH)
    bld["BA"] = norm_series(bld["BA_raw"], method=NORM_METHOD, p_low=P_LOW, p_high=P_HIGH)
    bld["SD"] = alpha_D * bld["KD"] + (1 - alpha_D) * bld["BA"]
    print(f"  SD computed — mean={bld['SD'].mean():.3f}, std={bld['SD'].std():.3f}")

    # ------------------------------------------------------------------
    # 5. SS: Distance Proxy (DP) + Orientation Proxy (OP)
    # ------------------------------------------------------------------

    def compute_SS_terms(i: int):
        """
        Returns (DP_raw_i, OP_raw_i, neighbor_count_i).
        DP_raw_i = mean(1 / (d_ij + eps))
        OP_raw_i = mean(g(theta_ij) / (d_ij + eps))
        where d_ij is wall-to-wall distance and g is the orientation weight.
        """
        Pi    = polys[i]
        phi_i = float(bld.iloc[i]["phi_deg"])
        idxs  = tree_query_indices(tree_polys, Pi.buffer(r_S))

        inv_sum = inv_orient_sum = 0.0
        m = 0
        for j in idxs:
            if j == i:
                continue
            dij = Pi.distance(polys[j])
            if dij > r_S:
                continue
            inv    = 1.0 / (dij + epsilon)
            phi_j  = float(bld.iloc[j]["phi_deg"])
            theta  = angle_difference_deg(phi_i, phi_j)
            orient = orientation_factor(theta, sigma_theta)
            inv_sum        += inv
            inv_orient_sum += orient * inv
            m += 1

        if m == 0:
            return 0.0, 0.0, 0
        return inv_sum / m, inv_orient_sum / m, m

    print("Computing SS components...")
    dp_raw = np.zeros(len(bld), dtype=float)
    op_raw = np.zeros(len(bld), dtype=float)
    m_cnt  = np.zeros(len(bld), dtype=int)
    for i in tqdm(range(len(bld)), desc="  SS (distance + orientation)"):
        dp_raw[i], op_raw[i], m_cnt[i] = compute_SS_terms(i)

    bld["DP_raw"]       = dp_raw
    bld["OP_raw"]       = op_raw
    bld["SS_neighbors"] = m_cnt

    bld["DP"] = norm_series(bld["DP_raw"], method=NORM_METHOD, p_low=P_LOW, p_high=P_HIGH)
    bld["OP"] = norm_series(bld["OP_raw"], method=NORM_METHOD, p_low=P_LOW, p_high=P_HIGH)
    bld["SS"] = alpha_S * bld["DP"] + (1 - alpha_S) * bld["OP"]
    print(f"  SS computed — mean={bld['SS'].mean():.3f}, std={bld['SS'].std():.3f}")

    # ------------------------------------------------------------------
    # 6. SSDD (linear and geometric blends) + provenance fields
    # ------------------------------------------------------------------
    bld["SSDD"]      = beta * bld["SD"] + (1 - beta) * bld["SS"]
    bld["SSDD_geom"] = (bld["SD"].clip(1e-9, 1) ** beta) * (bld["SS"].clip(1e-9, 1) ** (1 - beta))
    print(f"  SSDD computed — mean={bld['SSDD'].mean():.3f}, std={bld['SSDD'].std():.3f}")

    # Store all run parameters on each feature for self-describing outputs
    bld["SD_r_m"]   = float(r_D)
    bld["SS_r_m"]   = float(r_S)
    bld["alpha_D"]  = float(alpha_D)
    bld["alpha_S"]  = float(alpha_S)
    bld["beta"]     = float(beta)
    bld["sigma_th"] = float(sigma_theta)
    bld["eps_m"]    = float(epsilon)
    bld["kernel"]   = str(kernel_type)
    bld["w_area"]   = bool(weight_by_area)
    bld["norm"]     = str(NORM_METHOD)
    bld["p_low"]    = int(P_LOW)
    bld["p_high"]   = int(P_HIGH)

    # ------------------------------------------------------------------
    # 7. QA plots
    # ------------------------------------------------------------------
    if MAKE_QA_PLOTS:
        for col, title in [
            ("KD_raw", "KD_raw distribution"),
            ("BA_raw", "BA_raw (basal fraction) distribution"),
            ("SD",     "SD distribution"),
            ("DP_raw", "DP_raw distribution"),
            ("OP_raw", "OP_raw distribution"),
            ("SS",     "SS distribution"),
            ("SSDD",   "SSDD distribution"),
        ]:
            qa_hist(bld[col], title)

    # ------------------------------------------------------------------
    # 8. Save outputs
    # ------------------------------------------------------------------
    run_tag = (
        f"SDr{int(r_D)}_SSr{int(r_S)}_"
        f"aD{alpha_D:.2f}_aS{alpha_S:.2f}_b{beta:.2f}_"
        f"sig{sigma_theta:.0f}_eps{epsilon:.1f}_"
        f"norm{NORM_METHOD}"
    )
    out_gpkg = os.path.join(OUT_DIR, f"{RUN_NAME}_SSDD_{run_tag}.gpkg")
    out_csv  = os.path.join(OUT_DIR, f"{RUN_NAME}_SSDD_{run_tag}.csv")
    out_txt  = os.path.join(OUT_DIR, f"{RUN_NAME}_SSDD_{run_tag}_RunSummary.txt")

    out_gdf = bld.drop(columns=["rep_pt"], errors="ignore").copy()

    # GeoPackage (preferred; preserves geometry and field names)
    out_gdf.to_file(out_gpkg, layer="buildings_ssdd", driver="GPKG")

    # CSV (attribute table only)
    csv_cols = [
        "bld_id",
        "KD_raw", "BA_raw", "KD", "BA", "SD",
        "DP_raw", "OP_raw", "DP", "OP", "SS", "SS_neighbors",
        "SSDD", "SSDD_geom",
        "SD_r_m", "SS_r_m",
        "alpha_D", "alpha_S", "beta",
        "sigma_th", "eps_m",
        "kernel", "w_area",
        "norm", "p_low", "p_high",
    ]
    csv_cols_existing = [c for c in csv_cols if c in out_gdf.columns]
    out_gdf[csv_cols_existing].to_csv(out_csv, index=False)

    # Run summary text
    elapsed_sec = time.time() - t0
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("SSDD run summary\n")
        f.write("----------------\n")
        f.write(f"RUN_NAME: {RUN_NAME}\n")
        f.write(f"INPUT: {INPUT_PATH} (layer={INPUT_LAYER})\n")
        f.write(f"N buildings: {len(out_gdf):,}\n")
        f.write(f"CRS (analysis): {out_gdf.crs}\n\n")
        f.write("Parameters:\n")
        f.write(f"  SD radius r_D (m)    : {r_D}\n")
        f.write(f"  SS radius r_S (m)    : {r_S}\n")
        f.write(f"  alpha_D              : {alpha_D}\n")
        f.write(f"  alpha_S              : {alpha_S}\n")
        f.write(f"  beta                 : {beta}\n")
        f.write(f"  sigma_theta (deg)    : {sigma_theta}\n")
        f.write(f"  epsilon (m)          : {epsilon}\n")
        f.write(f"  kernel               : {kernel_type}\n")
        f.write(f"  weight_by_area       : {weight_by_area}\n")
        f.write(f"  normalization        : {NORM_METHOD} "
                f"(P_LOW={P_LOW}, P_HIGH={P_HIGH})\n\n")
        f.write("Outputs:\n")
        f.write(f"  GPKG: {out_gpkg} (layer=buildings_ssdd)\n")
        f.write(f"  CSV : {out_csv}\n\n")
        f.write(f"Elapsed seconds: {elapsed_sec:.2f}\n")

    print(f"\nWrote: {out_gpkg}")
    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_txt}")
    print(f"\nDone in {elapsed_sec:.1f}s")


if __name__ == "__main__":
    main()
