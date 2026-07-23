# %% [markdown]
# # 03 · Comparisons
#
# Includes the **4→2 collapse check**: does replacing the four production metrics
# {KD, BA, DP, OP} with two {SD, SS} lose univariate separation? Each old metric
# is a branch of a new one (KD/BA = SD unit/area weight; DP/OP = SS flat/oriented),
# so SD's grid ⊇ {KD, BA} and SS's grid ⊇ {DP, OP} by construction. The check
# reports, per concept, the best achievable separation in its branch alongside the
# real production metric — see `compare_collapse.csv`.
#
# **Primary** — two separation measures per form, per fire and pooled:
#
#   struct_auc   structure-level AUC = P(metric ranks a destroyed structure above
#                a survived one).  Rank-based, base-rate free.
#
#   nbhd_pauc    neighborhood pseudo-AUC = P(the form's block-mean ranks two
#                randomly chosen neighborhoods in the same order as their damage
#                rate).  Block-pair concordance (Somers' D mapped to [0,1]); 0.5 =
#                chance, same scale as AUC.  Reported alongside the Spearman
#                correlation between block-mean and block damage-rate.
#
# **Secondary** — fire profiles: where each form's variation lives.  Per fire, the
# structure-level distribution (median, IQR, CV) and the between/within-block
# variance split (ICC = share of variance that is between neighborhoods).
#
# Blocks are 1000 m UTM cells (≥ 5× the largest r_D, so between-block variance
# isn't an artifact of kernel overlap); neighborhood measures use blocks with
# ≥ 20 buildings.  Emits `compare_primary.csv` and `compare_profiles.csv`.

# %%
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


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
# Spatial blocking is chosen PER FIRE: a cell size from this ladder (all ≥3× the
# largest r_D so between-block variance isn't a kernel-overlap artifact). Sparse
# fires need coarser cells to hold enough buildings per neighborhood.
BLOCK_LADDER = (600, 1000, 1500, 2000, 3000, 4000)
MIN_BLOCK_N = 20


# %%
# ── load metrics + attach per-fire adaptive block id ─────────────────────────
metrics = pd.read_parquet(OUT_DIR / "metrics.parquet")
forms = pd.read_csv(OUT_DIR / "forms.csv", keep_default_na=False)
focal = pd.read_parquet(OUT_DIR / "focal.parquet")

metrics = metrics.merge(focal[["fire", "focal_id", "cent_x", "cent_y"]],
                        on=["fire", "focal_id"], how="left")
metric_cols = forms["col"].tolist()
FIRES = tuple(sorted(metrics["fire"].unique()))   # discovered from data — scales to new fires

# per fire, pick the cell size that yields the MOST usable (≥MIN_BLOCK_N) blocks —
# fine for dense fires, coarse for sparse ones — and namespace the id by fire.
def _bid(cx, cy, bs):
    return (cx // bs).astype(int).astype(str) + "_" + (cy // bs).astype(int).astype(str)

metrics["block"] = ""
block_choice = {}
for fire in FIRES:
    m = metrics["fire"] == fire
    cx, cy = metrics.loc[m, "cent_x"], metrics.loc[m, "cent_y"]
    best_bs, best_n = BLOCK_LADDER[0], -1
    for bs in BLOCK_LADDER:
        n = int((_bid(cx, cy, bs).value_counts() >= MIN_BLOCK_N).sum())
        if n > best_n:
            best_n, best_bs = n, bs
    block_choice[fire] = (best_bs, best_n)
    metrics.loc[m, "block"] = fire + "_" + _bid(cx, cy, best_bs)
print(f"{len(metrics):,} focals × {len(metric_cols)} forms · fires: {FIRES}")
print("adaptive block (size m, usable blocks):",
      {f: block_choice[f] for f in FIRES})


# %%
# ── measures ──────────────────────────────────────────────────────────────────
def auc(y: np.ndarray, s: np.ndarray) -> float:
    """AUC from midranks (== Mann-Whitney), handles ties correctly."""
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(s)
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def concordance(x: np.ndarray, r: np.ndarray) -> float:
    """Neighborhood pseudo-AUC: P(sign Δ block-mean == sign Δ damage-rate) + ½ ties.
    Block-pair concordance over pairs with different damage rate. [0,1], 0.5=chance."""
    x = np.asarray(x, float)
    r = np.asarray(r, float)
    if len(x) < 3:
        return float("nan")
    dx = np.sign(x[:, None] - x[None, :])
    dr = np.sign(r[:, None] - r[None, :])
    comp = dr != 0
    nc = comp.sum()
    if nc == 0:
        return float("nan")
    conc = ((dx == dr) & comp).sum()
    tie = ((dx == 0) & comp).sum()
    return (conc + 0.5 * tie) / nc


def var_decomp(s: np.ndarray, blocks: np.ndarray):
    """ANOVA split of a metric into between- vs within-block variance.
    Returns (total_CV, between_CV, within_CV, ICC) with ICC = between/total var."""
    gm = s.mean()
    if gm == 0 or not np.isfinite(gm):
        return (np.nan,) * 4
    df = pd.DataFrame({"s": s, "b": blocks})
    g = df.groupby("b")["s"]
    n_b = g.size().to_numpy()
    mean_b = g.mean().to_numpy()
    var_b = g.var(ddof=0).fillna(0.0).to_numpy()
    N = n_b.sum()
    between_var = (n_b * (mean_b - gm) ** 2).sum() / N
    within_var = (n_b * var_b).sum() / N
    total_var = between_var + within_var
    icc = between_var / total_var if total_var > 0 else np.nan
    cv = lambda v: np.sqrt(v) / gm
    return cv(total_var), cv(between_var), cv(within_var), icc


# %%
# ── primary + profiles, per fire ──────────────────────────────────────────────
prim_rows, prof_rows = [], []
meta = forms.set_index("col")[["family", "a", "b", "r"]].to_dict("index")

for fire in list(FIRES) + ["pooled"]:
    sub = metrics if fire == "pooled" else metrics[metrics["fire"] == fire]
    y = sub["destroyed"].to_numpy()
    blocks = sub["block"].to_numpy()

    # block-level tables computed once for all forms
    bstat = sub.groupby("block").agg(n=("destroyed", "size"),
                                     rate=("destroyed", "mean"))
    big = bstat.index[bstat["n"] >= MIN_BLOCK_N]
    bmean_all = sub.groupby("block")[metric_cols].mean()
    rate_big = bstat.loc[big, "rate"].to_numpy()

    for col in metric_cols:
        m = meta[col]
        s = sub[col].to_numpy()

        bm = bmean_all.loc[big, col].to_numpy()
        rho, _ = spearmanr(bm, rate_big) if len(big) >= 3 else (np.nan, None)
        prim_rows.append({
            "family": m["family"], "a": m["a"], "b": m["b"], "r": int(m["r"]),
            "form": col, "fire": fire,
            "struct_auc": auc(y, s),
            "nbhd_pauc": concordance(bm, rate_big),
            "spearman": float(rho),
            "n_blocks": int(len(big)),
        })

        tcv, bcv, wcv, icc = var_decomp(s, blocks)
        q25, med, q75 = np.quantile(s, [0.25, 0.50, 0.75])
        prof_rows.append({
            "family": m["family"], "a": m["a"], "b": m["b"], "r": int(m["r"]),
            "form": col, "fire": fire,
            "median": med, "q25": q25, "q75": q75,
            "total_cv": tcv, "between_cv": bcv, "within_cv": wcv, "icc_between": icc,
        })

primary = pd.DataFrame(prim_rows)
profiles = pd.DataFrame(prof_rows)
primary.to_csv(OUT_DIR / "compare_primary.csv", index=False)
profiles.to_csv(OUT_DIR / "compare_profiles.csv", index=False)
print(f"Wrote compare_primary.csv ({len(primary):,} rows) and "
      f"compare_profiles.csv ({len(profiles):,} rows)")


# %%
# ── quick look: pooled leaders per family ─────────────────────────────────────
pool = primary[primary["fire"] == "pooled"].copy()
for fam in ("SD", "SS"):
    top = (pool[pool["family"] == fam]
           .sort_values("struct_auc", ascending=False)
           .head(8)[["a", "b", "r", "struct_auc", "nbhd_pauc", "spearman"]])
    print(f"\n{fam}: top 8 by pooled structure-AUC")
    with pd.option_context("display.float_format", "{:.3f}".format):
        print(top.to_string(index=False))


# %% [markdown]
# ## 4 → 2 collapse check
#
# For each concept, the best univariate separation achievable within its branch,
# plus the real production metric where one exists. SD nests KD (unit) and BA
# (area); SS nests DP (flat) and OP (oriented).

# %%
# concept -> (new family, row filter within that family's grid)
BRANCH = {
    "KD": ("SD", lambda d: d["b"] == "unit"),
    "BA": ("SD", lambda d: d["b"] == "area"),
    "SD": ("SD", lambda d: d["b"].notna()),
    "DP": ("SS", lambda d: d["a"] == "flat"),
    "OP": ("SS", lambda d: d["a"] != "flat"),
    "SS": ("SS", lambda d: d["a"].notna()),
}
# production incumbent lookup: family=='incumbent', a in {KD,BA,DP,OP}
inc = primary[primary["family"] == "incumbent"]
prod = {(r.a, r.fire): (r.struct_auc, r.nbhd_pauc)
        for r in inc.itertuples(index=False)}

col_rows = []
for concept, (fam, filt) in BRANCH.items():
    for fire in list(FIRES) + ["pooled"]:
        d = primary[(primary["family"] == fam) & (primary["fire"] == fire)]
        d = d[filt(d)]
        if d.empty:
            continue
        bs = d.loc[d["struct_auc"].idxmax()]
        bn = d.loc[d["nbhd_pauc"].idxmax()]
        ps, pn = prod.get((concept, fire), (np.nan, np.nan))
        col_rows.append({
            "concept": concept, "fire": fire,
            "struct_best": bs["struct_auc"], "struct_form": bs["form"],
            "nbhd_best": bn["nbhd_pauc"], "nbhd_form": bn["form"],
            "production_struct": ps, "production_nbhd": pn,
        })
collapse = pd.DataFrame(col_rows)
collapse.to_csv(OUT_DIR / "compare_collapse.csv", index=False)
print(f"\nWrote compare_collapse.csv ({len(collapse)} rows)")

show = collapse[collapse["fire"] == "pooled"].set_index("concept")
print("\n4→2 collapse (pooled): best-achievable separation per concept")
with pd.option_context("display.float_format", "{:.3f}".format):
    print(show[["struct_best", "nbhd_best",
                "production_struct", "production_nbhd"]].to_string())
print("\nClaim: SD.best ≥ max(KD.best, BA.best);  SS.best ≥ max(DP.best, OP.best)")
for scale in ("struct_best", "nbhd_best"):
    sd = show.loc["SD", scale]; kd = show.loc["KD", scale]; ba = show.loc["BA", scale]
    ss = show.loc["SS", scale]; dp = show.loc["DP", scale]; op = show.loc["OP", scale]
    print(f"  {scale}:  SD {sd:.3f} vs KD {kd:.3f}/BA {ba:.3f}   |   "
          f"SS {ss:.3f} vs DP {dp:.3f}/OP {op:.3f}")
