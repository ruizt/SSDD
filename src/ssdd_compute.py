#!/usr/bin/env python3
"""CLI for computing raw SSDD metrics.

Reads a building footprint layer, optionally spatial-joins DINS structure
points, computes ``SD`` (structure density) and ``SS`` (structure
separation), and writes a tabular CSV plus a geometry-bearing GeoPackage keyed
on ``ssdd_id``.

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
        description="Compute raw SSDD metrics (SD, SS) for a building-footprint layer.",
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

    p.add_argument("--r-d", type=float, default=200.0,
                   help="SD neighborhood radius (m). From the metric-specification "
                        "experiments; 200 is the top of the swept grid.")
    p.add_argument("--r-s", type=float, default=50.0,
                   help="SS nearest-neighbour search radius (m). The selected form "
                        "is r_S-invariant over 25-75 m.")
    p.add_argument("--epsilon", type=float, default=0.5,
                   help="SS distance floor (m); bounds SS at 1/epsilon for touching "
                        "walls. Held fixed in the specification experiments.")
    p.add_argument("--kernel", default="uniform",
                   choices=["uniform", "epanechnikov", "quartic", "triweight"],
                   help="SD kernel shape. Immaterial to separation; uniform is the "
                        "simplest.")
    p.add_argument("--weight", default="root_area",
                   choices=["unit", "area", "root_area"],
                   help="SD neighbour weighting: 1, area, or sqrt(area).")
    p.add_argument("--agg", default="nn",
                   choices=["nn", "uniform", "power1", "power2"],
                   help="SS aggregation over neighbours. nn (default) is the "
                        "specified form; power1 reproduces the legacy DP/OP mean.")
    p.add_argument("--orient", default="flat",
                   choices=["flat", "gauss", "cos2", "cos4"],
                   help="SS orientation weight g(theta). flat (default) applies "
                        "no orientation preference.")
    p.add_argument("--sigma", type=float, default=10.0,
                   help="Gaussian orientation tolerance (deg); only used when "
                        "--orient gauss.")
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
        kernel=args.kernel,
        weight=args.weight,
        agg=args.agg,
        orient=args.orient,
        sigma=args.sigma,
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
        "ssdd_id", "bld_area", "phi_deg", "cent_x", "cent_y",
        "SD", "SS",
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
        f.write(f"  kernel={args.kernel}  weight={args.weight}\n")
        f.write(f"  agg={args.agg}  orient={args.orient}  "
                f"sigma={args.sigma}\n\n")
        f.write(f"Outputs:\n  {csv_path}\n  {gpkg_path}\n\n")
        f.write(f"Elapsed seconds: {elapsed:.2f}\n")

    print(f"\nWrote {csv_path}")
    print(f"Wrote {gpkg_path}")
    print(f"Wrote {log_path}")
    print(f"Done in {elapsed:.1f}s")


if __name__ == "__main__":
    main()
