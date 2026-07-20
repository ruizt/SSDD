"""Pure geometric helpers used by the SSDD metrics.

No file I/O, no spatial indexes — just math on shapely geometries / floats.
"""

from __future__ import annotations

import math
import warnings

KERNELS = ("uniform", "epanechnikov", "quartic", "triweight")


def kernel_value(u: float, kernel: str = "uniform") -> float:
    """Symmetric beta-family kernel K_n(u) = (1 - u^2)^n on u in [0, 1], else 0.

    n = 0 (uniform), 1 (Epanechnikov), 2 (quartic), 3 (triweight). ``uniform``
    is the default: the metric-specification experiments found kernel shape
    immaterial to separation (all forms within ~0.02 AUC and the per-fire
    winner varies), so the simplest — a plain truncated count — is preferred.
    """
    if u < 0.0 or u > 1.0:
        return 0.0
    if kernel == "uniform":
        return 1.0
    if kernel == "epanechnikov":
        return 1.0 - u * u
    if kernel == "quartic":
        return (1.0 - u * u) ** 2
    if kernel == "triweight":
        return (1.0 - u * u) ** 3
    raise ValueError(f"Unknown kernel type: {kernel!r}; expected one of {KERNELS}")


def dominant_orientation_degrees(poly) -> float:
    """Angle (deg, in [0, 180)) of the longest edge of the minimum rotated rectangle.

    The minimum rotated rectangle is the smallest-area enclosing rectangle at
    any angle; its longest edge points along the building's long axis, so this
    is the footprint's dominant orientation (undirected — 10° ≡ 190°).

    Not used by the shipped metrics (SS applies no orientation weighting), but
    kept as a per-building attribute for downstream analyses.

    Suppresses a benign ``RuntimeWarning: divide by zero encountered in
    oriented_envelope`` that shapely 2.1 emits for perfectly axis-aligned
    rectangles. The returned MRR — and therefore this function's output — is
    still numerically correct in that case.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                module=r"shapely\..*")
        mrr = poly.minimum_rotated_rectangle
    coords = list(mrr.exterior.coords)
    edges = []
    for k in range(4):
        x1, y1 = coords[k]
        x2, y2 = coords[k + 1]
        dx, dy = x2 - x1, y2 - y1
        edges.append((math.hypot(dx, dy), dx, dy))
    _, dx, dy = max(edges, key=lambda t: t[0])
    return math.degrees(math.atan2(dy, dx)) % 180.0


def angle_difference_deg(a: float, b: float) -> float:
    """Absolute orientation difference in degrees, folded into [0, 90]."""
    diff = abs(a - b) % 180.0
    diff = min(diff, 180.0 - diff)
    return min(diff, 90.0)
