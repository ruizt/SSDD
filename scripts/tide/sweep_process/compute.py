#!/usr/bin/env python3
"""Per-job entrypoint for the SSDD radius sweep on Tide.

Reads parameters from environment variables, runs ``compute_raw_metrics`` for
one ``(fire, r_D, r_S, r_NN)`` combination, spatial-joins DINS, and writes a
single CSV to ``$SSDD_OUT_DIR/<fire>_rD<r_D>_rS<r_S>/``.

Environment variables
---------------------
SSDD_FIRE        'eaton' or 'palisades'                                 (required)
SSDD_R_D         SD radius (m)                                          (required)
SSDD_R_S         SS radius (m)                                          (required)
SSDD_R_NN        Nearest-neighbor search radius (m)                     (default: 200)
SSDD_EPSILON     Distance floor (m)                                     (default: 0.5)
SSDD_SIGMA_THETA Orientation tolerance (deg)                            (default: 15)
SSDD_EPSG        Target CRS                                             (default: 32611)
SSDD_DATA_DIR    Root of input data (buildings/ and dins/ subdirs)      (default: /data)
SSDD_OUT_DIR     Output root; per-run subdir created underneath         (default: /jobs/output)

Test locally:
    SSDD_FIRE=palisades SSDD_R_D=100 SSDD_R_S=50 \\
      SSDD_DATA_DIR=$(pwd)/_data/raw SSDD_OUT_DIR=$(pwd)/_tmp \\
      python scripts/tide/sweep_process/compute.py
"""

from __future__ import annotations

import os
import time
from pathlib import Path

from ssdd.io import (
    ensure_projected_meters,
    join_dins,
    read_buildings,
    read_dins,
)
from ssdd.pipeline import RawMetricParams, compute_raw_metrics


def main() -> None:
    fire = os.environ["SSDD_FIRE"]
    r_D = float(os.environ["SSDD_R_D"])
    r_S = float(os.environ["SSDD_R_S"])
    r_NN = float(os.environ.get("SSDD_R_NN", "200"))
    epsilon = float(os.environ.get("SSDD_EPSILON", "0.5"))
    sigma_theta = float(os.environ.get("SSDD_SIGMA_THETA", "15"))
    epsg = int(os.environ.get("SSDD_EPSG", "32611"))
    data_dir = Path(os.environ.get("SSDD_DATA_DIR", "/data"))
    out_dir = Path(os.environ.get("SSDD_OUT_DIR", "/jobs/output"))

    fire_cap = fire[:1].upper() + fire[1:]
    buildings_path = data_dir / "buildings" / f"LARIAC6_Buildings_2020_{fire}.shp"
    dins_path = data_dir / "dins" / f"DINS_2025_{fire_cap}_Public_View.geojson"

    run_name = f"{fire}_rD{int(r_D)}_rS{int(r_S)}"
    run_dir = out_dir / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print(f"[ssdd-sweep] fire={fire}  r_D={r_D}  r_S={r_S}  r_NN={r_NN}", flush=True)
    print(f"[ssdd-sweep] reading buildings: {buildings_path}", flush=True)
    bld = read_buildings(str(buildings_path))
    bld = ensure_projected_meters(bld, epsg)
    print(f"[ssdd-sweep] N={len(bld):,} buildings; CRS={bld.crs}", flush=True)

    params = RawMetricParams(
        r_D=r_D, r_S=r_S, r_NN=r_NN,
        epsilon=epsilon, sigma_theta=sigma_theta,
    )
    print("[ssdd-sweep] computing raw metrics ...", flush=True)
    bld = compute_raw_metrics(bld, params=params, progress=False)

    print(f"[ssdd-sweep] reading DINS: {dins_path}", flush=True)
    dins = read_dins(str(dins_path))
    dins = ensure_projected_meters(dins, epsg)
    bld = join_dins(bld, dins, how="left")
    print(f"[ssdd-sweep] after DINS join: N={len(bld):,}", flush=True)

    csv_cols = [
        "ssdd_id", "bld_area", "phi_deg", "cent_x", "cent_y",
        "KD_raw", "BA_raw", "DP_raw", "OP_raw", "SS_neighbors",
        "dist_to_nearest_building", "bearing_to_nearest_building",
    ]
    csv_cols += [c for c in bld.columns if c not in csv_cols + ["geometry"]]

    csv_path = run_dir / f"{run_name}_raw_metrics.csv"
    bld.drop(columns="geometry").to_csv(csv_path, index=False, columns=csv_cols)

    print(f"[ssdd-sweep] wrote {csv_path}", flush=True)
    print(f"[ssdd-sweep] elapsed {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
