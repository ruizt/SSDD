# %% [markdown]
# # 02 · Structure-level metrics  (SD and SS)
#
# Two metric families, each a small grid of design choices:
#
# **SD — structure density** (collapses the old KD + BA). Kernel-weighted count of
# neighbors within `r_D`, over a kernel × weight grid:
#
#   SD = Σ_j  w_j · K(d_cent_j / r_D)  /  (π r_D²)
#
#   kernel K   uniform · Epanechnikov · quartic · triweight · exponential
#   weight w   unit (→ smoothed count)  ·  area (→ smoothed basal area)  ·  √area
#
# **SS — structure separation** (collapses the old DP + OP). Orientation-weighted
# proximity to neighbors within `r_S`, over an orientation × aggregation grid:
#
#   orientation g(θ)   flat (g≡1) · gaussian σ∈{5,10,20} · cos⁴θ · Lambert cos²θ
#   aggregation        nn      = g(θ*)·/(d*+ε)   at the nearest neighbor
#                      uniform = mean_j g(θ_j)                    (no distance weight)
#                      power1  = mean_j g(θ_j)/(d_j+ε)
#                      power2  = mean_j g(θ_j)/(d_j+ε)²
#
# With g≡1, SS nests the old DP forms exactly: flat·nn = DP nearest-neighbor,
# flat·power1 = incumbent DP. Both families are swept over radius grids.
#
# Emits `metrics.parquet` (wide: one row per focal, one column per form) and
# `forms.csv` (column → parsed family / kernel|orient / weight|agg / r).

# %%
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


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

OUT_DIR = ROOT / "dev/py/metric-specification/_out"

# ── design grids ──────────────────────────────────────────────────────────────
R_D_GRID = (50, 100, 150, 200)
R_S_GRID = (25, 50, 75)
EPS = 0.5
EXP_H = 0.25            # exponential KD decay length (in units of u = d/r_D)

SD_KERNELS = {
    "uniform":      lambda u: np.ones_like(u),            # n=0
    "sqrt":         lambda u: np.sqrt(1.0 - u * u),       # n=1/2  (slower than Epa)
    "epanechnikov": lambda u: 1.0 - u * u,                # n=1
    "quartic":      lambda u: (1.0 - u * u) ** 2,         # n=2
    "triweight":    lambda u: (1.0 - u * u) ** 3,         # n=3
    "exponential":  lambda u: np.exp(-u / EXP_H),         # faster than triweight
}
SD_WEIGHTS = {
    "unit":      lambda a: np.ones_like(a),
    "area":      lambda a: a,
    "root_area": lambda a: np.sqrt(a),
}
SS_ORIENT = {
    "flat":    lambda th: np.ones_like(th),   # g≡1: no orientation preference
    "gauss5":  lambda th: np.exp(-(th / 5.0) ** 2),
    "gauss10": lambda th: np.exp(-(th / 10.0) ** 2),
    "gauss20": lambda th: np.exp(-(th / 20.0) ** 2),
    "cos4":    lambda th: np.cos(np.radians(th)) ** 4,
    "lambert": lambda th: np.cos(np.radians(th)) ** 2,
}
SS_AGG = ("nn", "uniform", "power1", "power2")


# %%
# ── load primitives, factorize focal groups ──────────────────────────────────
focal = pd.read_parquet(OUT_DIR / "focal.parquet")
pairs = pd.read_parquet(OUT_DIR / "pairs.parquet")
print(f"{len(focal):,} focals · {len(pairs):,} pairs")

key = pairs["fire"].astype(str) + "|" + pairs["focal_id"].astype(str)
codes, uniques = pd.factorize(key)
NG = len(uniques)

d_cent = pairs["d_cent"].to_numpy(float)
d_wall = pairs["d_wall"].to_numpy(float)
phi    = pairs["phi_diff"].to_numpy(float)      # already folded to [0, 90]
n_area = pairs["n_area"].to_numpy(float)
invd1  = 1.0 / (d_wall + EPS)
invd2  = invd1 * invd1


def gsum(x: np.ndarray) -> np.ndarray:
    return np.bincount(codes, weights=x, minlength=NG)


cols: dict[str, np.ndarray] = {}
meta: list[dict] = []


# %%
# ── SD: kernel × weight × r_D ─────────────────────────────────────────────────
for r_D in R_D_GRID:
    in_rd = (d_cent <= r_D).astype(float)
    u = np.clip(d_cent / r_D, 0.0, 1.0)
    norm = math.pi * r_D * r_D
    for kname, kfun in SD_KERNELS.items():
        K = kfun(u) * in_rd
        for wname, wfun in SD_WEIGHTS.items():
            s = gsum(K * wfun(n_area)) / norm
            col = f"SD__{kname}|{wname}__r{r_D}"
            cols[col] = s
            meta.append({"col": col, "family": "SD", "a": kname, "b": wname, "r": r_D})
print(f"SD forms: {sum(m['family']=='SD' for m in meta)}")


# %%
# ── SS: orientation × aggregation × r_S ───────────────────────────────────────
for r_S in R_S_GRID:
    in_rs = (d_wall <= r_S).astype(float)
    cnt = gsum(in_rs)
    have = cnt > 0

    # nearest neighbor within r_S: one pair-row per focal group
    dm = np.where(in_rs > 0, d_wall, np.inf)
    order = np.lexsort((dm, codes))
    sc = codes[order]
    first = np.empty(len(order), bool)
    first[0] = True
    first[1:] = sc[1:] != sc[:-1]
    nn_rows = order[first]                      # candidate nearest row per code
    nn_code = codes[nn_rows]
    nn_ok = np.isfinite(dm[nn_rows])            # nearest actually within r_S

    for oname, ofun in SS_ORIENT.items():
        g = ofun(phi)
        num_u = gsum(g * in_rs)
        num_1 = gsum(g * invd1 * in_rs)
        num_2 = gsum(g * invd2 * in_rs)
        for agg in SS_AGG:
            v = np.zeros(NG)
            if agg == "uniform":
                np.divide(num_u, cnt, out=v, where=have)
            elif agg == "power1":
                np.divide(num_1, cnt, out=v, where=have)
            elif agg == "power2":
                np.divide(num_2, cnt, out=v, where=have)
            elif agg == "nn":
                contrib = (g[nn_rows] * invd1[nn_rows])
                v[nn_code[nn_ok]] = contrib[nn_ok]
            col = f"SS__{oname}|{agg}__r{r_S}"
            cols[col] = v
            meta.append({"col": col, "family": "SS", "a": oname, "b": agg, "r": r_S})
print(f"SS forms: {sum(m['family']=='SS' for m in meta)}")


# %%
# ── assemble wide table over ALL focals (fill 0 for focals absent from pairs) ──
out = pd.concat([pd.DataFrame({"key": uniques}), pd.DataFrame(cols)], axis=1)

# incumbent production metrics (from the gpkgs, carried in focal.parquet) — the
# real 4-metric baseline for the collapse check. Production radii: KD/BA r_D=100,
# DP/OP r_S=50.
INCUMBENT = {"KD_raw": ("incumbent", "KD", 100),
             "BA_raw": ("incumbent", "BA", 100),
             "DP_raw": ("incumbent", "DP", 50),
             "OP_raw": ("incumbent", "OP", 50)}
inc_cols = [c for c in INCUMBENT if c in focal.columns]

base = focal[["fire", "focal_id", "destroyed", *inc_cols]].copy()
base["key"] = base["fire"].astype(str) + "|" + base["focal_id"].astype(str)
metrics = base.merge(out, on="key", how="left").drop(columns="key")
metric_cols = list(cols.keys())
metrics[metric_cols] = metrics[metric_cols].fillna(0.0)

for c in inc_cols:
    fam, a, r = INCUMBENT[c]
    meta.append({"col": c, "family": fam, "a": a, "b": "production", "r": r})
metric_cols = metric_cols + inc_cols

metrics.to_parquet(OUT_DIR / "metrics.parquet", index=False)
forms = pd.DataFrame(meta)
forms.to_csv(OUT_DIR / "forms.csv", index=False)
print(f"Wrote metrics.parquet ({metrics.shape[0]:,} × {len(metric_cols)} forms) "
      f"and forms.csv")


# %%
# ── nesting sanity checks: new forms must reproduce the old KD/DP numbers ──────
def auc(y, s):
    pos, neg = s[y == 1], s[y == 0]
    u, _ = mannwhitneyu(pos, neg, alternative="greater")
    return u / (len(pos) * len(neg))


y = metrics["destroyed"].to_numpy()
checks = {
    "SD__uniform|unit__r200  (== old KD uniform, 0.691)": "SD__uniform|unit__r200",
    "SS__flat|power1__r50    (== old DP incumbent, 0.536)": "SS__flat|power1__r50",
    "SS__flat|nn__r50        (== old nn_proxy, 0.610)": "SS__flat|nn__r50",
    "SD__quartic|unit__r200  (== old KD quartic, 0.678)": "SD__quartic|unit__r200",
}
print("\nnesting checks (pooled AUC):")
for label, col in checks.items():
    print(f"  {label:<52s} -> {auc(y, metrics[col].to_numpy()):.3f}")
