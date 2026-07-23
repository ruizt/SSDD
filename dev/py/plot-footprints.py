# %% [markdown]
# # Footprint figures
#
# Two products, both driven off the fires present in `focal.parquet` (so they
# scale automatically as fires are added):
#
#   block_map_{fire}.png   per-fire footprints with the 600 m analysis grid
#                          overlaid; blocks that enter the neighborhood measures
#                          (≥20 buildings) are shaded.
#   footprints_bare.png    all fires side by side, DINS-labeled footprints only,
#                          shared vertical scale, transparent background — for
#                          dropping onto a slide.
#
# Reads the processed gpkgs directly; independent of the metric-spec pipeline
# and of the ssdd package (fires discovered from `_data/processed/`).

# %%
from __future__ import annotations

import warnings
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle

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
except NameError:
    ROOT = _find_root(Path.cwd())

OUT_DIR = ROOT / "dev/py/_out"          # footprint figures live outside metric-specification
OUT_DIR.mkdir(parents=True, exist_ok=True)
BLOCK = 600.0            # must match 03_compare.BLOCK
MIN_BLOCK_N = 20         # must match 03_compare.MIN_BLOCK_N

CONTEXT   = "#d9dde2"
SURVIVED  = "#2b7a8c"
DESTROYED = "#c94b34"
GRID      = "#7a828c"
USABLE    = "#f2d98a"
BARE      = "#3f4a5c"

# fires discovered from the processed gpkgs — no dependency on the metric pipeline
FIRES = sorted(p.parent.name for p in ROOT.glob("_data/processed/*/*_buildings.gpkg"))
print(f"fires: {FIRES}")


def load_fire(fire: str) -> gpd.GeoDataFrame:
    g = gpd.read_file(ROOT / f"_data/processed/{fire}/{fire}_buildings.gpkg")
    if g.crs is None or g.crs.to_epsg() != 32611:
        g = g.to_crs(32611)
    return g


# %% [markdown]
# ## Per-fire block-overlay maps

# %%
def block_map(fire: str) -> None:
    g = load_fire(fire)
    lab = g[g["DAMAGE"].notna()].copy()
    lab["destroyed"] = (lab["DAMAGE"] == "Destroyed (>50%)").astype(int)
    cent = lab.geometry.centroid
    lab["bx"] = (cent.x // BLOCK).astype(int)
    lab["by"] = (cent.y // BLOCK).astype(int)

    blk = lab.groupby(["bx", "by"]).agg(n=("destroyed", "size"),
                                        n_d=("destroyed", "sum")).reset_index()
    blk["usable"] = blk["n"] >= MIN_BLOCK_N
    n_usable = int(blk["usable"].sum())

    lxmin, lymin, lxmax, lymax = lab.total_bounds
    _, fymin, _, fymax = g.total_bounds
    M = 250.0
    xmin, xmax = lxmin - M, lxmax + M
    ymin, ymax = fymin - 100.0, fymax + 100.0
    aspect = (xmax - xmin) / (ymax - ymin)
    fig, ax = plt.subplots(figsize=(8.0 * aspect, 8.0))

    for _, r in blk[blk["usable"]].iterrows():
        ax.add_patch(Rectangle((r["bx"] * BLOCK, r["by"] * BLOCK), BLOCK, BLOCK,
                               facecolor=USABLE, edgecolor="none", alpha=0.45, zorder=0))
    g.plot(ax=ax, facecolor=CONTEXT, edgecolor="none", zorder=1)
    lab[lab["destroyed"] == 0].plot(ax=ax, facecolor=SURVIVED, edgecolor="none", zorder=2)
    lab[lab["destroyed"] == 1].plot(ax=ax, facecolor=DESTROYED, edgecolor="none", zorder=3)

    x0 = np.floor(xmin / BLOCK) * BLOCK
    y0 = np.floor(ymin / BLOCK) * BLOCK
    for x in np.arange(x0, xmax + BLOCK, BLOCK):
        ax.plot([x, x], [y0, ymax], color=GRID, lw=0.6, alpha=0.6, zorder=4)
    for y in np.arange(y0, ymax + BLOCK, BLOCK):
        ax.plot([x0, xmax], [y, y], color=GRID, lw=0.6, alpha=0.6, zorder=4)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.set_xlabel("UTM 11N easting (m)")
    ax.set_ylabel("UTM 11N northing (m)")
    n_d = int(lab["destroyed"].sum())
    ax.set_title(f"{fire.capitalize()} — footprints & {BLOCK:.0f} m analysis grid\n"
                 f"{len(g):,} footprints · {len(lab):,} DINS-labeled "
                 f"({n_d:,} destroyed / {len(lab) - n_d:,} survived) · "
                 f"{n_usable} usable blocks (≥{MIN_BLOCK_N})", fontsize=11)
    ax.legend(handles=[
        Patch(facecolor=DESTROYED, label="destroyed (DINS)"),
        Patch(facecolor=SURVIVED, label="survived (DINS)"),
        Patch(facecolor=CONTEXT, label="unlabeled footprint"),
        Patch(facecolor=USABLE, alpha=0.45, label=f"usable block (≥{MIN_BLOCK_N})"),
    ], loc="upper right", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    png = OUT_DIR / f"block_map_{fire}.png"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png.name}  ({n_usable} usable blocks)")


for fire in FIRES:
    block_map(fire)


# %% [markdown]
# ## Bare transparent footprints (all fires, shared scale)

# %%
labeled, dims = {}, {}
for fire in FIRES:
    lab = load_fire(fire)
    lab = lab[lab["DAMAGE"].notna()].copy()
    labeled[fire] = lab
    xmin, ymin, xmax, ymax = lab.total_bounds
    dims[fire] = (xmax - xmin, ymax - ymin)

ratios = [dims[f][0] / dims[f][1] for f in FIRES]
fig_h = 4.0
fig, axes = plt.subplots(1, len(FIRES),
                         figsize=(fig_h * sum(ratios) * 1.05, fig_h),
                         gridspec_kw={"width_ratios": ratios, "wspace": 0.04},
                         squeeze=False)
for ax, fire in zip(axes[0], FIRES):
    labeled[fire].plot(ax=ax, facecolor=BARE, edgecolor="none")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.margins(0)
png = OUT_DIR / "footprints_bare.png"
fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.02, transparent=True)
plt.close(fig)
print(f"wrote {png.name}")
