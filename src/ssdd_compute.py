#!/usr/bin/env python3
"""CLI for computing raw SSDD metrics.

Reads a building footprint layer, optionally spatial-joins DINS structure
points, computes ``KD_raw``, ``BA_raw``, ``DP_raw``, ``OP_raw``, and writes a
tabular CSV plus a geometry-bearing GeoPackage keyed on ``ssdd_id``.

The package deliberately stops at raw metrics — normalization, blending and
predictive modeling are downstream choices.

Run from ``SSDD/src/``::

    python ssdd_compute.py \\
        --buildings ../_data/raw/buildings/LARIAC6_Buildings_2020_eaton.shp \\
        --dins      ../_data/raw/dins/DINS_2025_Eaton_Public_View.geojson \\
        --output    ../_data/processed/eaton \\
        --run-name  eaton

Pass ``--dins-only`` to drop buildings without a DINS hit (the LA-fires
"burned subset" filter).
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from ssdd.io import (
    ensure_projected_meters,
    join_dins,
    read_buildings,
    read_dins,
)
from ssdd.pipeline import RawMetricParams, compute_raw_metrics


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute raw SSDD metrics for a building-footprint layer.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--buildings", required=True, help="Path to building footprints (SHP/GPKG/GeoJSON).")
    p.add_argument("--buildings-layer", default=None, help="Layer name (GPKG/FGDB).")
    p.add_argument("--dins", default=None, help="Optional DINS point file to spatial-join.")
    p.add_argument("--dins-layer", default=None, help="DINS layer name (GPKG/FGDB).")
    p.add_argument("--dins-only", action="store_true",
                   help="Drop buildings with no DINS match (inner join).")
    p.add_argument("--output", required=True, help="Output directory.")
    p.add_argument("--run-name", default="run", help="Prefix for output files.")
    p.add_argument("--epsg", type=int, default=32611,
                   help="Target CRS (projected, meters). Default: UTM 11N.")

    p.add_argument("--r-d", type=float, default=100.0, help="SD buffer radius (m).")
    p.add_argument("--r-s", type=float, default=50.0, help="SS search radius (m).")
    p.add_argument("--epsilon", type=float, default=0.5, help="Distance floor (m).")
    p.add_argument("--sigma-theta", type=float, default=15.0, help="Orientation tolerance (deg).")
    p.add_argument("--kernel", default="quartic", help="Kernel name.")
    p.add_argument("--weight-by-area", action="store_true",
                   help="Weight KD neighbor contributions by footprint area.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print(f"Reading buildings: {args.buildings}")
    bld = read_buildings(args.buildings, layer=args.buildings_layer)
    print(f"  Input CRS: {bld.crs}  N={len(bld):,}")
    bld = ensure_projected_meters(bld, args.epsg)
    print(f"  Analysis CRS: {bld.crs}")

    params = RawMetricParams(
        r_D=args.r_d,
        r_S=args.r_s,
        epsilon=args.epsilon,
        sigma_theta=args.sigma_theta,
        kernel=args.kernel,
        weight_by_area=args.weight_by_area,
    )

    print("Computing raw metrics...")
    bld = compute_raw_metrics(bld, params=params)

    if args.dins:
        print(f"Reading DINS: {args.dins}")
        dins = read_dins(args.dins, layer=args.dins_layer)
        dins = ensure_projected_meters(dins, args.epsg)
        print(f"  DINS points: {len(dins):,}")
        how = "inner" if args.dins_only else "left"
        bld = join_dins(bld, dins, how=how)
        print(f"  After DINS join ({how}): N={len(bld):,}")

    csv_path = out_dir / f"{args.run_name}_raw_metrics.csv"
    gpkg_path = out_dir / f"{args.run_name}_buildings.gpkg"
    log_path = out_dir / f"{args.run_name}_compute_log.txt"

    csv_cols = [
        "ssdd_id", "bld_area", "phi_deg",
        "KD_raw", "BA_raw", "DP_raw", "OP_raw", "SS_neighbors",
    ]
    csv_cols += [c for c in bld.columns if c not in csv_cols + ["geometry"]]
    bld.drop(columns="geometry").to_csv(csv_path, index=False, columns=csv_cols)
    bld.to_file(gpkg_path, layer="buildings_raw", driver="GPKG")

    elapsed = time.time() - t0
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("SSDD Stage 1 — compute log\n")
        f.write("----------------------------\n")
        f.write(f"run_name        : {args.run_name}\n")
        f.write(f"buildings       : {args.buildings}\n")
        f.write(f"dins            : {args.dins}\n")
        f.write(f"dins_only       : {args.dins_only}\n")
        f.write(f"n buildings out : {len(bld):,}\n")
        f.write(f"analysis CRS    : EPSG:{args.epsg}\n\n")
        f.write("Parameters:\n")
        f.write(f"  r_D={args.r_d}  r_S={args.r_s}  epsilon={args.epsilon}\n")
        f.write(f"  sigma_theta={args.sigma_theta}  kernel={args.kernel}\n")
        f.write(f"  weight_by_area={args.weight_by_area}\n\n")
        f.write(f"Outputs:\n  {csv_path}\n  {gpkg_path}\n\n")
        f.write(f"Elapsed seconds: {elapsed:.2f}\n")

    print(f"\nWrote {csv_path}")
    print(f"Wrote {gpkg_path}")
    print(f"Wrote {log_path}")
    print(f"Done in {elapsed:.1f}s")


if __name__ == "__main__":
    main()
