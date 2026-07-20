# %% [markdown]
# # 04 · Plots
#
# All figures for the SD/SS sensitivity study, driven by the CSVs from 03 (plus
# the raw pairs for the orientation histogram). Each cell writes one PNG to
# `_out/` (transparent, 300 dpi). Runs top-to-bottom or cell-by-cell in Positron.
#
#   sd_kernels.png      SD kernel shapes (beta n=0..3 + exponential)
#   ss_orient.png       SS orientation weights (flat, gaussians, cos⁴, Lambert)
#   ss_distance.png     SS distance aggregations (nn / uniform / 1/dᵖ)
#   sd_heat_*.png       SD kernel×weight AUC heatmaps, faceted by r_D
#   ss_heat_*.png       SS orient×agg AUC heatmaps, faceted by r_S
#   scale_scatter.png   structure-AUC vs neighborhood-pseudo-AUC, all forms
#   profiles_icc.png    between-block variance share (ICC) for selected forms

# %%
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


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

OUT = ROOT / "dev/py/metric-specification/_out"

primary = pd.read_csv(OUT / "compare_primary.csv", keep_default_na=False)
profiles = pd.read_csv(OUT / "compare_profiles.csv", keep_default_na=False)
collapse = pd.read_csv(OUT / "compare_collapse.csv", keep_default_na=False)
for df in (primary, profiles, collapse):
    for c in ("struct_auc", "nbhd_pauc", "spearman", "r",
              "median", "q25", "q75", "total_cv", "between_cv",
              "within_cv", "icc_between",
              "struct_best", "nbhd_best", "production_struct", "production_nbhd"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

SD_KERN = ["uniform", "sqrt", "epanechnikov", "quartic", "triweight", "exponential"]
SD_WEIGHT = ["unit", "area", "root_area"]
SS_ORIENT = ["flat", "gauss5", "gauss10", "gauss20", "cos4", "lambert"]
SS_AGG = ["nn", "uniform", "power1", "power2"]
R_D_GRID = [50, 100, 150, 200]
R_S_GRID = [25, 50, 75]
SD_R_TABLE = 200          # radius slice shown in the per-fire forms table
SS_R_TABLE = 50

# fires discovered from the data (excludes the pooled aggregate) — scales as
# more fires are added upstream; nothing here is hard-coded to two fires.
FIRES = [f for f in sorted(primary["fire"].unique()) if f != "pooled"]
FIRE_COLS = FIRES + ["pooled"]
FIRE_COLORS = {f: c for f, c in zip(FIRES, plt.cm.Dark2.colors)}

INK = "#2d3748"


def save(fig, name):
    p = OUT / name
    fig.savefig(p, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {p.name}")


# %% [markdown]
# ## SD kernel shapes

# %%
u = np.linspace(0, 1.25, 700)
EXP_H = 0.25


def beta(u, n):
    return np.where(u <= 1, np.maximum(0.0, 1 - np.clip(u, 0, 1) ** 2) ** n, 0.0)


kd_forms = [
    ("$n=0$  uniform",         beta(u, 0),   "#9fb3c8", "-"),
    ("$n=\\frac{1}{2}$  sqrt", beta(u, 0.5), "#5b8def", "-"),
    ("$n=1$  Epanechnikov",    beta(u, 1),   "#2f6f9f", "-"),
    ("$n=2$  quartic",         beta(u, 2),   "#274b6d", "-"),
    ("$n=3$  triweight",       beta(u, 3),   "#16263a", "-"),
    ("exponential  $h=0.25$",  np.where(u <= 1, np.exp(-u / EXP_H), 0.0),
     "#c94b34", (0, (1.4, 1.6))),
]

fig, ax = plt.subplots(figsize=(7.2, 3.8))
for label, k, c, ls in kd_forms:
    ax.plot(u, k, lw=2.2, color=c, ls=ls, label=label)
ax.axvline(1.0, color="#3f4a5c", lw=0.9, ls=(0, (2, 3)), alpha=0.7)
ax.annotate("$u=1$  ($d=r_D$)", xy=(1.0, 1.02), xytext=(-3, 0),
            textcoords="offset points", fontsize=8, color="#3f4a5c", ha="right", va="top")
ax.set_xlim(0, 1.25); ax.set_ylim(0, 1.06)
ax.set_xlabel("$u = d\\,/\\,r_D$"); ax.set_ylabel("$K(u)$")
ax.set_title("SD kernel shapes  ($K_n(u)=(1-u^2)^n$ + exponential)", fontsize=11, pad=8)
ax.grid(True, ls=":", lw=0.5, alpha=0.45)
ax.legend(fontsize=8.5, loc="center left", bbox_to_anchor=(1.01, 0.5), borderaxespad=0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
save(fig, "sd_kernels.png")


# %% [markdown]
# ## SS orientation weights

# %%
theta = np.linspace(0, 90, 900)
pairs = pd.read_parquet(OUT / "pairs.parquet", columns=["d_wall", "phi_diff"])
dtheta = pairs.loc[pairs["d_wall"] <= 50, "phi_diff"].to_numpy()

ss_orient_forms = [
    ("flat  $g\\equiv 1$",            np.ones_like(theta),                 "#94a3b8", (0, (1, 2))),
    ("gaussian  $\\sigma=5°$",        np.exp(-(theta / 5.0) ** 2),         "#8c2f1c", "-"),
    ("gaussian  $\\sigma=10°$",       np.exp(-(theta / 10.0) ** 2),        "#c94b34", "-"),
    ("gaussian  $\\sigma=20°$",       np.exp(-(theta / 20.0) ** 2),        "#e8926b", "-"),
    ("$\\cos^4\\theta$",              np.cos(np.radians(theta)) ** 4,      "#6c4f9c", (0, (3, 2))),
    ("Lambert  $\\cos^2\\theta$",     np.cos(np.radians(theta)) ** 2,      "#2b7a8c", (0, (5, 2.5))),
]

fig, ax = plt.subplots(figsize=(7.0, 3.8))
ax2 = ax.twinx()
ax2.hist(dtheta, bins=np.arange(0, 92, 2), density=True, color="#3f4a5c", alpha=0.12)
ax2.set_yticks([]); ax2.spines[:].set_visible(False)
for label, g, c, ls in ss_orient_forms:
    ax.plot(theta, g, lw=2.1, color=c, ls=ls, label=label)
ax.set_xlim(0, 90); ax.set_ylim(0, 1.05); ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
ax.set_xlabel("orientation difference $\\Delta\\theta$ (deg)"); ax.set_ylabel("$g(\\Delta\\theta)$")
ax.set_title("SS orientation weights", fontsize=11, pad=8)
ax.grid(True, ls=":", lw=0.5, alpha=0.45)
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
ax.annotate("parallel", xy=(2, 0.05), fontsize=8, color="#3f4a5c", ha="left")
ax.annotate("perpendicular", xy=(88, 0.05), fontsize=8, color="#3f4a5c", ha="right")
ax.legend(fontsize=8.5, loc="center left", bbox_to_anchor=(1.08, 0.5), borderaxespad=0)
for s in ("top",):
    ax.spines[s].set_visible(False)
fig.tight_layout()
save(fig, "ss_orient.png")


# %% [markdown]
# ## SS distance aggregations

# %%
EPS = 0.5
d = np.linspace(0, 50, 1200)
PAIR_Q25, PAIR_MED, PAIR_Q75 = 18.2, 30.9, 40.5
DMIN = 1.27


def normd(g, at):
    return g / at


ss_dist_forms = [
    ("power1  $1/(d{+}\\epsilon)$",  normd(1 / (d + EPS), 1 / (DMIN + EPS)),   "#c94b34", "-"),
    ("power2  $1/(d{+}\\epsilon)^2$", normd(1 / (d + EPS) ** 2, 1 / (DMIN + EPS) ** 2), "#e08a5d", "-"),
]
fig, ax = plt.subplots(figsize=(6.6, 3.8))
ax.axvspan(PAIR_Q25, PAIR_Q75, color="#3f4a5c", alpha=0.10)
ax.axvline(PAIR_MED, color="#3f4a5c", lw=1.0, ls=(0, (2, 3)), alpha=0.8)
ax.annotate("median neighbour", xy=(PAIR_MED, 4e-4), xytext=(-6, 0),
            textcoords="offset points", fontsize=8, color="#3f4a5c", va="bottom", ha="right")
for label, g, c, ls in ss_dist_forms:
    ax.plot(d, g, lw=2.2, color=c, ls=ls, label=label)
# uniform (flat within r_S) + nn markers, schematic
ax.axhline(1.0, color="#94a3b8", lw=0.7, alpha=0.6)
ax.plot([], [], color="#2b7a8c", lw=2.0, ls=(0, (4, 2)), label="uniform  (flat within $r_S$)")
ax.plot([], [], "o", color=INK, ms=6, label="nn  (nearest neighbour only)")
ax.set_yscale("log"); ax.set_xlim(0, 50); ax.set_ylim(1e-4, 20)
ax.set_xlabel("wall-to-wall distance $d$ (m)"); ax.set_ylabel("weight $/\\,$weight$(d_{min})$")
ax.set_title("SS distance aggregations", fontsize=11, pad=8)
ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.4)
ax.legend(fontsize=8.5, loc="upper right", title=f"$\\epsilon={EPS:g}$ m", title_fontsize=8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
save(fig, "ss_distance.png")


# %% [markdown]
# ## Heatmaps — AUC over the design grid (pooled)

# %%
def heat_facets(family, rowvals, colvals, rgrid, measure, fname, title):
    pool = primary[(primary["fire"] == "pooled") & (primary["family"] == family)]
    dmin, dmax = pool[measure].min(), pool[measure].max()
    # sequential when the family never runs backwards; diverging when it crosses 0.5
    if dmin >= 0.5:
        cmap = "YlGn"
        norm = mcolors.Normalize(vmin=dmin - 0.003, vmax=dmax)
    else:
        cmap = "RdYlGn"
        norm = mcolors.TwoSlopeNorm(vcenter=0.5, vmin=min(dmin, 0.46), vmax=dmax)
    cmo = plt.get_cmap(cmap)

    def txt_color(v):
        r, g, b, _ = cmo(norm(v))
        return "white" if (0.299 * r + 0.587 * g + 0.114 * b) < 0.5 else INK

    n = len(rgrid)
    fig, axes = plt.subplots(1, n, figsize=(3.05 * n, 0.62 * len(rowvals) + 1.4),
                             squeeze=False)
    axes = axes[0]
    im = None
    best = pool.loc[pool[measure].idxmax()]
    for ax, rv in zip(axes, rgrid):
        sub = pool[pool["r"] == rv]
        M = (sub.pivot(index="a", columns="b", values=measure)
                .reindex(index=rowvals, columns=colvals))
        im = ax.imshow(M.values, cmap=cmap, norm=norm, aspect="auto")
        for i in range(len(rowvals)):
            for j in range(len(colvals)):
                v = M.values[i, j]
                if np.isnan(v):
                    continue
                is_best = (best["a"] == rowvals[i] and best["b"] == colvals[j]
                           and best["r"] == rv)
                ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                        fontsize=8, color=txt_color(v),
                        fontweight="bold" if is_best else "normal")
                if is_best:
                    ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                               edgecolor=INK, lw=2.2))
        ax.set_xticks(range(len(colvals)))
        ax.set_xticklabels(colvals, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(len(rowvals)))
        ax.set_yticklabels(rowvals if ax is axes[0] else [], fontsize=8)
        rlab = "r_D" if family == "SD" else "r_S"
        ax.set_title(f"${rlab}={rv}$", fontsize=9.5)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    fig.suptitle(title, fontsize=11.5, y=1.02)
    fig.tight_layout(rect=(0, 0, 0.94, 1))
    cax = fig.add_axes([0.955, 0.18, 0.012, 0.62])
    fig.colorbar(im, cax=cax).ax.tick_params(labelsize=7.5)
    save(fig, fname)


heat_facets("SD", SD_KERN, SD_WEIGHT, R_D_GRID, "struct_auc",
            "sd_heat_struct.png", "SD  ·  structure-level AUC   (kernel × weight)")
heat_facets("SD", SD_KERN, SD_WEIGHT, R_D_GRID, "nbhd_pauc",
            "sd_heat_nbhd.png", "SD  ·  neighborhood pseudo-AUC   (kernel × weight)")
heat_facets("SS", SS_ORIENT, SS_AGG, R_S_GRID, "struct_auc",
            "ss_heat_struct.png", "SS  ·  structure-level AUC   (orientation × aggregation)")
heat_facets("SS", SS_ORIENT, SS_AGG, R_S_GRID, "nbhd_pauc",
            "ss_heat_nbhd.png", "SS  ·  neighborhood pseudo-AUC   (orientation × aggregation)")


# %% [markdown]
# ## Two-scale summary — structure AUC vs neighborhood pseudo-AUC

# %%
pool = primary[primary["fire"] == "pooled"].copy()
fig, ax = plt.subplots(figsize=(5.6, 5.2))
colors = {"SD": "#2f6f9f", "SS": "#c94b34"}
for fam, c in colors.items():
    s = pool[pool["family"] == fam]
    ax.scatter(s["struct_auc"], s["nbhd_pauc"], s=26, color=c, alpha=0.55,
               edgecolor="white", lw=0.4, label=fam)
lo, hi = 0.44, 0.76
ax.plot([lo, hi], [lo, hi], color="#94a3b8", lw=0.8, ls="--", zorder=0)
ax.axvline(0.5, color="#cbd5e1", lw=0.6, zorder=0)
ax.axhline(0.5, color="#cbd5e1", lw=0.6, zorder=0)
# label the leaders
for _, row in pool.sort_values("struct_auc").tail(1).iterrows():
    ax.annotate("SD uniform·√area·200", (row["struct_auc"], row["nbhd_pauc"]),
                fontsize=8, color=colors["SD"], xytext=(-4, 6),
                textcoords="offset points", ha="right")
nn = pool[(pool["family"] == "SS") & (pool["b"] == "nn") & (pool["r"] == 50)].iloc[0]
ax.annotate("SS flat·nn", (nn["struct_auc"], nn["nbhd_pauc"]),
            fontsize=8, color=colors["SS"], xytext=(6, -2), textcoords="offset points")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.set_xlabel("structure-level AUC"); ax.set_ylabel("neighborhood pseudo-AUC")
ax.set_title("Where the signal lives (pooled)", fontsize=11, pad=8)
ax.legend(fontsize=9, loc="lower right", title="family", title_fontsize=9)
ax.grid(True, ls=":", lw=0.5, alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
save(fig, "scale_scatter.png")


# %% [markdown]
# ## Fire profiles — between-block variance share (ICC) for selected forms

# %%
SELECT = [
    ("SD uniform·√area·200", "SD", "uniform", "root_area", 200),
    ("SD uniform·unit·200",  "SD", "uniform", "unit", 200),
    ("SD quartic·area·100",  "SD", "quartic", "area", 100),
    ("SS flat·nn",           "SS", "flat", "nn", 50),
    ("SS gauss5·power2·25",  "SS", "gauss5", "power2", 25),
]
labels = [s[0] for s in SELECT]
x = np.arange(len(SELECT))
nf = len(FIRES)
w = 0.8 / nf
fig, ax = plt.subplots(figsize=(7.6, 4.0))
for k, fire in enumerate(FIRES):
    icc = []
    for _, fam, a, b, r in SELECT:
        row = profiles[(profiles["fire"] == fire) & (profiles["family"] == fam)
                       & (profiles["a"] == a) & (profiles["b"] == b) & (profiles["r"] == r)]
        icc.append(row["icc_between"].iloc[0] if len(row) else np.nan)
    ax.bar(x + (k - (nf - 1) / 2) * w, icc, w, color=FIRE_COLORS[fire],
           label=fire, edgecolor="white")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8.5)
ax.set_ylabel("ICC = between-block variance share")
ax.set_title("Where each form's variation lives  (higher = more neighborhood-level)",
             fontsize=10.5, pad=8)
ax.set_ylim(0, 1)
ax.legend(fontsize=9, title="fire", title_fontsize=9)
ax.grid(True, axis="y", ls=":", lw=0.5, alpha=0.45)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
save(fig, "profiles_icc.png")


# %% [markdown]
# ## 4 → 2 collapse — best-achievable separation per concept

# %%
cc = collapse[collapse["fire"] == "pooled"].set_index("concept")
ORDER = ["KD", "BA", "SD", "DP", "OP", "SS"]
XPOS = {"KD": 0, "BA": 1, "SD": 2, "DP": 3.7, "OP": 4.7, "SS": 5.7}
FACE = {"KD": "#9dc3e0", "BA": "#9dc3e0", "SD": "#2f6f9f",
        "DP": "#eaa89a", "OP": "#eaa89a", "SS": "#c94b34"}
NEW = {"SD", "SS"}

fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4), sharey=True)
for ax, (bestcol, prodcol, ttl) in zip(
        axes,
        [("struct_best", "production_struct", "structure-level AUC"),
         ("nbhd_best", "production_nbhd", "neighborhood pseudo-AUC")]):
    for cpt in ORDER:
        x = XPOS[cpt]
        v = cc.loc[cpt, bestcol]
        ax.bar(x, v - 0.5, bottom=0.5, width=0.82, color=FACE[cpt],
               edgecolor=INK if cpt in NEW else "none",
               lw=2.0 if cpt in NEW else 0, zorder=2)
        ax.text(x, v + 0.004, f"{v:.3f}", ha="center", va="bottom",
                fontsize=8.5, fontweight="bold" if cpt in NEW else "normal", color=INK)
        p = cc.loc[cpt, prodcol]
        if np.isfinite(p):
            ax.plot(x, p, "o", ms=6, color=INK, zorder=3)
            ax.plot([x - 0.42, x + 0.42], [p, p], color=INK, lw=1.0, zorder=3)
    ax.axhline(0.5, color="#94a3b8", lw=1.0, zorder=1)
    ax.set_xticks([XPOS[c] for c in ORDER])
    ax.set_xticklabels(ORDER, fontsize=10)
    ax.set_ylim(0.42, 0.76)
    ax.set_title(ttl, fontsize=11, pad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    ax.text(1.0, 0.432, "density → SD", ha="center", fontsize=8.5, color="#2f6f9f")
    ax.text(4.7, 0.432, "separation → SS", ha="center", fontsize=8.5, color="#c94b34")

axes[0].set_ylabel("best achievable (bar) · production ●")
fig.suptitle("Collapsing 4 metrics into 2 loses no univariate separation",
             fontsize=12.5, y=1.02)
fig.tight_layout()
save(fig, "collapse.png")


# %% [markdown]
# ## Per-fire result tables (both comparisons)
# Plots above are pooled; these tables disaggregate by fire. Columns are
# discovered from the data, so adding fires just widens the tables.

# %%
HDR, BAD = "#2d3748", "#c94b34"


def draw_table(col_labels, row_labels, values, fname, title,
               bold_rows=None, box_best=True, fmt="{:.3f}", figw=None):
    bold_rows = bold_rows or set()
    values = np.asarray(values, float)
    nrow, ncol = values.shape
    text = [[row_labels[i]] + [fmt.format(values[i, j]) for j in range(ncol)]
            for i in range(nrow)]
    fig, ax = plt.subplots(figsize=(figw or (1.6 + 1.0 * ncol), 0.34 * nrow + 1.2))
    ax.axis("off")
    tbl = ax.table(cellText=text, colLabels=[""] + col_labels,
                   cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.32)
    for j in range(ncol + 1):
        c = tbl[0, j]
        c.set_facecolor(HDR)
        c.set_text_props(color="white", fontweight="bold")
    best = values.argmax(axis=0)
    for i in range(nrow):
        for j in range(ncol + 1):
            tbl[i + 1, j].set_edgecolor("#dee2e6")
            tbl[i + 1, j].set_facecolor("#eef1f5" if i in bold_rows else "white")
        tbl[i + 1, 0].set_text_props(
            ha="left", color=INK,
            fontweight="bold" if i in bold_rows else "normal")
        for j in range(ncol):
            props = {}
            if i in bold_rows or (box_best and best[j] == i):
                props["fontweight"] = "bold"
            if values[i, j] < 0.5:
                props["color"] = BAD
            if props:
                tbl[i + 1, j + 1].set_text_props(**props)
    ax.set_title(title, fontsize=11.5, pad=12, color=INK)
    fig.tight_layout()
    save(fig, fname)


# number of metrics — collapse, per fire
concepts = ["KD", "BA", "DP", "OP", "SD", "SS"]
cl = [f"{f}\nstruct" for f in FIRE_COLS] + [f"{f}\nnbhd" for f in FIRE_COLS]
V = np.full((len(concepts), 2 * len(FIRE_COLS)), np.nan)
for i, cpt in enumerate(concepts):
    for k, fire in enumerate(FIRE_COLS):
        r = collapse[(collapse["concept"] == cpt) & (collapse["fire"] == fire)]
        if len(r):
            V[i, k] = r["struct_best"].iloc[0]
            V[i, len(FIRE_COLS) + k] = r["nbhd_best"].iloc[0]
draw_table(cl, concepts, V, "table_num_metrics.png",
           "Number of metrics — best univariate separation by fire\n"
           "(SD/SS = the two shipped metrics; KD/BA/DP/OP = their branches)",
           bold_rows={4, 5}, box_best=False, figw=1.7 + 1.05 * 2 * len(FIRE_COLS))


# functional forms — per fire, at the pooled-winning radius
def forms_table(fam, kerns, weights, r_fix, fname, title):
    d = primary[(primary["family"] == fam) & (primary["r"] == r_fix)]
    rl, rows = [], []
    for a in kerns:
        for b in weights:
            rl.append(f"{a}·{b}")
            row = []
            for fire in FIRE_COLS:
                s = d[(d["a"] == a) & (d["b"] == b) & (d["fire"] == fire)]
                row.append(s["struct_auc"].iloc[0] if len(s) else np.nan)
            for fire in FIRE_COLS:
                s = d[(d["a"] == a) & (d["b"] == b) & (d["fire"] == fire)]
                row.append(s["nbhd_pauc"].iloc[0] if len(s) else np.nan)
            rows.append(row)
    cl = [f"{f}\nstruct" for f in FIRE_COLS] + [f"{f}\nnbhd" for f in FIRE_COLS]
    draw_table(cl, rl, np.array(rows), fname, title,
               box_best=True, figw=2.2 + 1.05 * 2 * len(FIRE_COLS))


forms_table("SD", SD_KERN, SD_WEIGHT, SD_R_TABLE, "table_forms_sd.png",
            f"Functional forms — SD (kernel·weight, r_D={SD_R_TABLE}) by fire\n"
            "bold = best in column;  full radius sweep in compare_primary.csv")
forms_table("SS", SS_ORIENT, SS_AGG, SS_R_TABLE, "table_forms_ss.png",
            f"Functional forms — SS (orientation·aggregation, r_S={SS_R_TABLE}) by fire\n"
            "bold = best in column;  full radius sweep in compare_primary.csv")


# %% [markdown]
# ## SD & SS value distributions per fire (destroyed vs survived)

# %%
mt = pd.read_parquet(OUT / "metrics.parquet")
DESTROYED, SURVIVED = "#c94b34", "#2b7a8c"
DIST_FORMS = [("SD__uniform|unit__r200", "SD  (count density, $r_D{=}200$)"),
              ("SS__flat|nn__r50",       "SS  (nearest-neighbour, $r_S{=}50$)")]

fig, axes = plt.subplots(len(DIST_FORMS), len(FIRES),
                         figsize=(3.2 * len(FIRES), 5.0), squeeze=False)
for r, (form, flabel) in enumerate(DIST_FORMS):
    pos_all = mt[form][mt[form] > 0]
    lo, hi = pos_all.quantile(0.01), pos_all.max()
    bins = np.logspace(np.log10(lo), np.log10(hi), 30)
    for c, fire in enumerate(FIRES):
        ax = axes[r][c]
        sub = mt[mt["fire"] == fire]
        v = sub[form].to_numpy(); y = sub["destroyed"].to_numpy()
        for outcome, color, lab in [(0, SURVIVED, "survived"), (1, DESTROYED, "destroyed")]:
            vv = v[y == outcome]
            ax.hist(vv[vv > 0], bins=bins, density=True, histtype="stepfilled",
                    color=color, alpha=0.42, lw=1.4, edgecolor=color, label=lab)
        zt = np.mean(v == 0) * 100
        ax.set_xscale("log"); ax.set_xlim(lo, hi)
        ax.set_yticks([])
        if zt >= 1:
            ax.text(0.02, 0.95, f"{zt:.0f}% zero", transform=ax.transAxes,
                    fontsize=7.5, va="top", color="#3f4a5c")
        if r == 0:
            ax.set_title(fire, fontsize=11)
        if c == 0:
            ax.set_ylabel(flabel, fontsize=9.5)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        if r == 0 and c == len(FIRES) - 1:
            ax.legend(fontsize=8, loc="upper right", frameon=False)
fig.suptitle("Metric distributions by fire — destroyed vs survived", fontsize=12.5, y=1.0)
fig.tight_layout()
save(fig, "sd_ss_distributions.png")

print("\nall figures written")
