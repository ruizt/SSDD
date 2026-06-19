"""SSDD — building-level raw metrics for the Structure Separation Distance Density work.

The package computes four raw per-building metrics from a footprint layer:

  Structure Density family (SD)
    KD_raw   kernel-weighted neighbor count per area, point-based.
    BA_raw   fraction of nearby ground covered by buildings, polygon-based.

  Structure Separation family (SS)
    DP_raw   mean inverse wall-to-wall distance to nearby neighbors.
    OP_raw   DP weighted by orientation alignment between focal and neighbor.

Normalization, blending and downstream modeling are intentionally out of scope —
they are the user's analysis to design. See ``dev/scripts/py/ssdd.py`` for the
historical convex-blending prototype.

Submodules
----------
io         File I/O, CRS handling, and the optional DINS spatial join.
geometry   Pure helpers (orientation, kernel, angle).
metrics    The four raw metric implementations.
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
    orientation_factor,
    quartic_kernel,
)
from .metrics import (
    compute_KD_series,
    compute_BA_series,
    compute_SS_terms_df,
)
from .pipeline import compute_raw_metrics
from . import synthetic

__all__ = [
    "read_buildings",
    "read_dins",
    "ensure_projected_meters",
    "join_dins",
    "dominant_orientation_degrees",
    "angle_difference_deg",
    "orientation_factor",
    "quartic_kernel",
    "compute_KD_series",
    "compute_BA_series",
    "compute_SS_terms_df",
    "compute_raw_metrics",
    "synthetic",
]
