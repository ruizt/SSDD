"""Pure geometric helpers used by the SSDD metrics.

No file I/O, no spatial indexes — just math on shapely geometries / floats.
"""

from __future__ import annotations

import math
import warnings


def quartic_kernel(u: float) -> float:
    """K(u) = (1 - u^2)^2 for u in [0, 1], else 0."""
    if u < 0.0 or u > 1.0:
        return 0.0
    return (1.0 - u * u) ** 2


def kernel_value(u: float, kernel: str = "quartic") -> float:
    if kernel == "quartic":
        return quartic_kernel(u)
    raise ValueError(f"Unknown kernel type: {kernel}")


def dominant_orientation_degrees(poly) -> float:
    """Angle (deg, in [0, 180)) of the longest edge of the minimum rotated rectangle.

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


def orientation_factor(theta_deg: float, sigma_theta: float) -> float:
    """Gaussian-style orientation weight g(theta) = exp(-(theta / sigma)^2)."""
    if sigma_theta <= 0:
        return 1.0
    return math.exp(-((theta_deg / sigma_theta) ** 2))
