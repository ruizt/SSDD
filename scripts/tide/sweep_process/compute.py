#!/usr/bin/env python3
"""Per-job entrypoint for the SSDD radius sweep on Tide.

Reads one fire's processed buildings gpkg (DINS damage already joined), and for
one ``r_D`` computes the **candidate metric forms** the CV sweep compares:

  SD (4)  {uniform, quartic} kernel × {root_area, unit} weight   (all depend on r_D)
  SS (3)  nearest-neighbour, orient ∈ {flat, gauss, cos2}        (r_D-independent)

r_S is gone: SS uses the true nearest neighbour with no isolation cutoff, so the
sweep's only metric-radius dimension is r_D. All forms are computed in a single
vectorised pass over shared STRtrees (bulk ``dwithin`` for SD, bulk
``query_nearest`` for SS).

Writes one CSV to ``$SSDD_OUT_DIR/<fire>_rD<r_D>/`` with ssdd_id, cent_x, cent_y,
DAMAGE, and the 7 metric columns.

Environment variables
---------------------
SSDD_FIRE      fire name (matches _data/processed/<fire>/<fire>_buildings.gpkg)  (required)
SSDD_R_D       SD radius (m)                                                     (required)
SSDD_EPSILON   SS distance floor (m)                                            (default 0.5)
SSDD_SIGMA     Gaussian orientation tolerance (deg), for SS_gauss               (default 10)
SSDD_EPSG      Target CRS                                                       (default 32611)
SSDD_CHUNK     Sources per SD neighbour-query batch (memory cap)                (default 4000)
SSDD_DATA_DIR  Root holding processed/<fire>/<fire>_buildings.gpkg              (default /data)
SSDD_OUT_DIR   Output root                                                      (default /jobs/output)

Test locally:
    SSDD_FIRE=mountain SSDD_R_D=200 \\
      SSDD_DATA_DIR=$(pwd)/_data/processed SSDD_OUT_DIR=$(pwd)/_tmp \\
      python scripts/tide/sweep_process/compute.py
"""

from __future__ import annotations

import math
import os
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.strtree import STRtree

from ssdd.geometry import dominant_orientation_degrees

# Sources per SD neighbour-query batch. Caps peak pair-array memory; see the
# SD block in main(). Lower it if a very dense fire still OOMs.
CHUNK = int(os.environ.get("SSDD_CHUNK", "4000"))

SD_KERNELS = ("uniform", "quartic")
SD_WEIGHTS = ("root_area", "unit")
SS_ORIENTS = ("flat", "gauss", "cos2")


def _kernel(u: np.ndarray, name: str) -> np.ndarray:
    if name == "uniform":
        return (u <= 1.0).astype(float)
    if name == "quartic":
        return np.where(u <= 1.0, (1.0 - u * u) ** 2, 0.0)
    raise ValueError(name)


def _fold(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    d = np.abs(a - b) % 180.0
    d = np.minimum(d, 180.0 - d)
    return np.minimum(d, 90.0)


def _orient(theta: np.ndarray, name: str, sigma: float) -> np.ndarray:
    if name == "flat":
        return np.ones_like(theta)
    if name == "gauss":
        return np.exp(-((theta / sigma) ** 2))
    if name == "cos2":
        return np.cos(np.radians(theta)) ** 2
    raise ValueError(name)


def main() -> None:
    fire = os.environ["SSDD_FIRE"]
    r_D = float(os.environ["SSDD_R_D"])
    eps = float(os.environ.get("SSDD_EPSILON", "0.5"))
    sigma = float(os.environ.get("SSDD_SIGMA", "10"))
    epsg = int(os.environ.get("SSDD_EPSG", "32611"))
    data_dir = Path(os.environ.get("SSDD_DATA_DIR", "/data"))
    out_dir = Path(os.environ.get("SSDD_OUT_DIR", "/jobs/output"))

    gpkg = data_dir / fire / f"{fire}_buildings.gpkg"
    run_name = f"{fire}_rD{int(r_D)}"
    run_dir = out_dir / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print(f"[ssdd-sweep] fire={fire} r_D={r_D}", flush=True)
    bld = gpd.read_file(gpkg)
    if bld.crs is None or bld.crs.to_epsg() != epsg:
        bld = bld.to_crs(epsg)
    bld = bld[bld.geometry.notna() & ~bld.geometry.is_empty].reset_index(drop=True)
    n = len(bld)
    print(f"[ssdd-sweep] N={n:,} buildings", flush=True)

    polys = bld.geometry.values
    cents = bld.geometry.centroid.values
    cx = np.array([c.x for c in cents]); cy = np.array([c.y for c in cents])
    areas = bld.geometry.area.to_numpy()
    phi = np.array([dominant_orientation_degrees(g) for g in polys])
    out = {"ssdd_id": bld.get("ssdd_id", np.arange(n)),
           "cent_x": cx, "cent_y": cy,
           "DAMAGE": bld.get("DAMAGE")}

    # ── SD: chunked dwithin passes, then 4 weightings ────────────────────────
    # The neighbour-pair array is the memory bottleneck: it grows as N·ρ·πr_D²,
    # so a dense fire at a large r_D (eaton at r_D=300) OOMs if queried in one
    # shot. Querying CHUNK sources at a time caps peak pair memory independent
    # of N and r_D; the bincount accumulates across chunks, so results are
    # identical to the single-pass version.
    norm = math.pi * r_D * r_D
    sd = {f"SD_{kn}_{wn}": np.zeros(n) for kn in SD_KERNELS for wn in SD_WEIGHTS}
    tc = STRtree(cents)
    for s in range(0, n, CHUNK):
        idx = np.arange(s, min(s + CHUNK, n))
        qi, ti = tc.query(cents[idx], predicate="dwithin", distance=r_D)
        gi = idx[qi]                       # query index is chunk-local -> globalise
        d = np.hypot(cx[gi] - cx[ti], cy[gi] - cy[ti])
        keep = (gi != ti) & (d <= r_D)
        gi, ti, d = gi[keep], ti[keep], d[keep]
        u = d / r_D
        aj = areas[ti]
        wmap = {"unit": np.ones_like(aj), "root_area": np.sqrt(aj)}
        for kn in SD_KERNELS:
            K = _kernel(u, kn)
            for wn in SD_WEIGHTS:
                sd[f"SD_{kn}_{wn}"] += np.bincount(gi, weights=wmap[wn] * K, minlength=n)
    for k, v in sd.items():
        out[k] = v / norm

    # ── SS: one bulk nearest pass, then 3 orientation weights ────────────────
    tp = STRtree(polys)
    (bi, nj), nd = tp.query_nearest(polys, exclusive=True, all_matches=False,
                                    return_distance=True)
    theta = _fold(phi[bi], phi[nj])
    for orient in SS_ORIENTS:
        g = _orient(theta, orient, sigma)
        col = np.zeros(n)
        col[bi] = g / (nd + eps)        # buildings absent from bi are alone -> 0
        out[f"SS_{orient}"] = col

    import pandas as pd
    csv_path = run_dir / f"{run_name}_metrics.csv"
    pd.DataFrame(out).to_csv(csv_path, index=False)
    print(f"[ssdd-sweep] wrote {csv_path}  ({time.time() - t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
