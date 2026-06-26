#!/usr/bin/env python3
"""collect.py — assemble per-job sweep CSVs into one long-format table.

Run after fetch.sh has pulled per-job outputs to ``RAW_DIR``. Each per-job
CSV is read, tagged with ``fire``, ``r_D``, ``r_S`` parsed from its
directory name, and concatenated into a single CSV at ``OUT_FILE``.

Usage (from the repo root, after fetch.sh):
    python scripts/tide/sweep_process/collect.py

Optional environment overrides:
    RAW_DIR   directory holding per-job subdirs (default: _data/processed/sweep)
    OUT_FILE  output CSV path (default: _data/processed/sweep/sweep_all.csv)
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd

RAW_DIR = Path(os.environ.get("RAW_DIR", "_data/processed/sweep"))
OUT_FILE = Path(os.environ.get("OUT_FILE", "_data/processed/sweep/sweep_all.csv"))

NAME_RE = re.compile(r"^(?P<fire>[a-z]+)_rD(?P<r_D>\d+)_rS(?P<r_S>\d+)$")


def collect() -> pd.DataFrame:
    if not RAW_DIR.exists():
        raise SystemExit(f"RAW_DIR not found: {RAW_DIR}")

    frames = []
    for sub in sorted(RAW_DIR.iterdir()):
        if not sub.is_dir():
            continue
        m = NAME_RE.match(sub.name)
        if not m:
            print(f"  skipping {sub.name} (doesn't match <fire>_rD<n>_rS<n>)")
            continue
        csv = sub / f"{sub.name}_raw_metrics.csv"
        if not csv.exists():
            print(f"  WARNING: {csv} missing")
            continue
        df = pd.read_csv(csv)
        df.insert(0, "fire", m["fire"])
        df.insert(1, "r_D", int(m["r_D"]))
        df.insert(2, "r_S", int(m["r_S"]))
        frames.append(df)
        print(f"  + {sub.name}: {len(df):,} rows")

    if not frames:
        raise SystemExit(f"No matching per-job CSVs found under {RAW_DIR}")

    return pd.concat(frames, ignore_index=True)


def main() -> None:
    print(f"Collecting sweep outputs from: {RAW_DIR}")
    out = collect()
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_FILE, index=False)
    print(f"\nWrote {OUT_FILE}  ({len(out):,} rows, {len(out.columns)} cols)")
    print("Combinations present:")
    print(
        out.groupby(["fire", "r_D", "r_S"])
        .size()
        .rename("n_rows")
        .reset_index()
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
