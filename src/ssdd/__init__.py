"""SSDD — building-level raw metrics for the Structure Separation & Density work.

The package computes two raw per-building metrics from a footprint layer, as
selected by the metric-specification experiments
(``dev/py/metric-specification/``):

  SD   structure density — kernel-weighted neighbour mass per unit area
           within ``r_D`` (uniform kernel, root-area weighting by default).
           Replaces the earlier KD + BA pair: the weighting choice subsumes
           the basal-area idea.

  SS   structure separation — inverse wall-to-wall distance to the nearest
           neighbour within ``r_S`` (no averaging, no orientation weighting).
           Replaces the earlier DP + OP pair: nearest-neighbour distance
           dominates every averaged / orientation-weighted alternative.

Both are coded so higher = more vulnerable (denser / closer). Normalization,
blending and downstream modeling are intentionally out of scope — they are the
user's analysis to design.

Submodules
----------
io         File I/O, CRS handling, and the optional DINS spatial join.
geometry   Pure helpers (kernel, orientation, angle folding).
metrics    The two raw metric implementations.
pipeline   End-to-end orchestration: ``compute_raw_metrics``.
synthetic  Parametric geometry generators for tests and sensitivity work.
"""

from .io import (
    read_buildings,
    read_dins,
    ensure_projected_meters,
    join_dins,
)
from .geometry import (
    dominant_orientation_degrees,
    angle_difference_deg,
    kernel_value,
)
from .metrics import (
    compute_SD_series,
    compute_SS_series,
)
from .pipeline import RawMetricParams, compute_raw_metrics
from . import synthetic

__all__ = [
    "read_buildings",
    "read_dins",
    "ensure_projected_meters",
    "join_dins",
    "dominant_orientation_degrees",
    "angle_difference_deg",
    "kernel_value",
    "compute_SD_series",
    "compute_SS_series",
    "RawMetricParams",
    "compute_raw_metrics",
    "synthetic",
]
