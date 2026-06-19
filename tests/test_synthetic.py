"""Accuracy tests against analytically known synthetic geometries.

Each test builds a small fixture from :mod:`ssdd.synthetic`, runs
:func:`ssdd.pipeline.compute_raw_metrics`, and checks the raw metric columns
against hand-computed expected values.

The math behind each expectation is spelled out in the docstring of the test
so it doubles as documentation of how the metrics behave on simple geometry.
"""

from __future__ import annotations

import math

import pytest

from ssdd import synthetic
from ssdd.geometry import orientation_factor
from ssdd.pipeline import RawMetricParams, compute_raw_metrics


def _buffer_window_area(poly, r_D: float) -> float:
    """Area of buffer(poly, r_D) as shapely actually computes it.

    The analytical Minkowski-sum area ``A + perimeter * r_D + pi * r_D^2``
    is *very close* but not identical, because shapely's default ``buffer``
    approximates the circle with 16 segments per quadrant. We test against
    the implementation shapely uses, since that's what the pipeline calls.
    """
    return poly.buffer(r_D).area


# ----------------------------------------------------------------------------
# Isolated building
# ----------------------------------------------------------------------------

def test_isolated_building_has_zero_neighbors():
    """One 10x20 rectangle, no neighbors.

    Expected:
      KD_raw = 0           (no other buildings)
      DP_raw = 0           (no neighbors)
      OP_raw = 0           (no neighbors)
      SS_neighbors = 0
      BA_raw = area / Minkowski-window-area, self-only.
    """
    gdf = synthetic.isolated_building(width=10.0, length=20.0)
    params = RawMetricParams(r_D=100.0, r_S=50.0, epsilon=0.5, sigma_theta=15.0)
    out = compute_raw_metrics(gdf, params=params, progress=False)

    assert out["KD_raw"].iloc[0] == 0.0
    assert out["DP_raw"].iloc[0] == 0.0
    assert out["OP_raw"].iloc[0] == 0.0
    assert int(out["SS_neighbors"].iloc[0]) == 0
    assert out["bld_area"].iloc[0] == pytest.approx(200.0)

    win_area = _buffer_window_area(out.geometry.iloc[0], params.r_D)
    expected_ba = 200.0 / win_area
    assert out["BA_raw"].iloc[0] == pytest.approx(expected_ba, rel=1e-9)


# ----------------------------------------------------------------------------
# Parallel pair
# ----------------------------------------------------------------------------

def test_parallel_pair_recovers_inverse_distance():
    """Two 10x20 rectangles, walls 10 m apart, both axis-aligned.

    Each building has exactly one neighbor at wall distance 10 m, so:
      SS_neighbors = 1
      DP_raw = 1 / (10 + epsilon)
      Both share orientation -> g(0) = 1 -> OP_raw = DP_raw
      Center-to-center distance = spacing + width = 20 m
      KD_raw = K(20/100) / (pi * 100^2), with K_quartic(0.2) = (1 - 0.04)^2 = 0.9216
    """
    spacing = 10.0
    width, length = 10.0, 20.0
    gdf = synthetic.pair(spacing=spacing, orientation_offset_deg=0.0,
                         width=width, length=length)
    params = RawMetricParams(r_D=100.0, r_S=50.0, epsilon=0.5, sigma_theta=15.0)
    out = compute_raw_metrics(gdf, params=params, progress=False)

    assert (out["SS_neighbors"] == 1).all()

    expected_dp = 1.0 / (spacing + params.epsilon)
    assert out["DP_raw"].iloc[0] == pytest.approx(expected_dp, rel=1e-9)
    assert out["DP_raw"].iloc[1] == pytest.approx(expected_dp, rel=1e-9)

    # Parallel => OP == DP
    assert out["OP_raw"].iloc[0] == pytest.approx(expected_dp, rel=1e-9)
    assert out["OP_raw"].iloc[1] == pytest.approx(expected_dp, rel=1e-9)

    u = (spacing + width) / params.r_D
    expected_kd = (1.0 - u * u) ** 2 / (math.pi * params.r_D ** 2)
    assert out["KD_raw"].iloc[0] == pytest.approx(expected_kd, rel=1e-3)


# ----------------------------------------------------------------------------
# Perpendicular pair
# ----------------------------------------------------------------------------

def test_perpendicular_pair_attenuates_OP():
    """Two rectangles, B2 rotated 90 deg. Same DP_raw, but OP_raw suppressed.

    Expected:
      DP_raw same shape as the parallel test, using the *actual* wall-to-wall
      distance after rotation (which is < spacing because the rotated corner
      reaches inward — see synthetic.pair docstring).
      OP_raw = g(theta) * DP_raw where theta = angle_difference_deg(phi_i, phi_j)
      folded into [0, 90]. With sigma_theta=15, g(90) = exp(-(90/15)^2) ~= 2.3e-16.
    """
    spacing = 30.0  # big enough that rotation doesn't bring walls into contact
    width, length = 10.0, 20.0
    gdf = synthetic.pair(spacing=spacing, orientation_offset_deg=90.0,
                         width=width, length=length)
    params = RawMetricParams(r_D=100.0, r_S=50.0, epsilon=0.5, sigma_theta=15.0)
    out = compute_raw_metrics(gdf, params=params, progress=False)

    # Use the actual geometry to get true wall-to-wall distance.
    actual_dist = gdf.geometry.iloc[0].distance(gdf.geometry.iloc[1])
    expected_dp = 1.0 / (actual_dist + params.epsilon)
    assert out["DP_raw"].iloc[0] == pytest.approx(expected_dp, rel=1e-9)

    expected_g = orientation_factor(90.0, params.sigma_theta)
    expected_op = expected_g * expected_dp
    assert out["OP_raw"].iloc[0] == pytest.approx(expected_op, abs=1e-12)
    # 90-degree case should drive OP essentially to zero
    assert out["OP_raw"].iloc[0] < 1e-10


# ----------------------------------------------------------------------------
# Touching buildings — saturation of inverse-distance
# ----------------------------------------------------------------------------

def test_touching_pair_saturates_at_one_over_epsilon():
    """Two rectangles sharing a wall.

    Wall-to-wall distance is 0, so:
      DP_raw = 1 / (0 + epsilon) = 1 / epsilon
      OP_raw = g(0) * DP_raw = DP_raw (parallel by default)
    """
    eps = 0.5
    gdf = synthetic.touching_pair(orientation_offset_deg=0.0)
    params = RawMetricParams(r_D=100.0, r_S=50.0, epsilon=eps, sigma_theta=15.0)
    out = compute_raw_metrics(gdf, params=params, progress=False)

    assert out["DP_raw"].iloc[0] == pytest.approx(1.0 / eps, rel=1e-9)
    assert out["OP_raw"].iloc[0] == pytest.approx(1.0 / eps, rel=1e-9)


# ----------------------------------------------------------------------------
# 3x3 grid — neighbor counts at known radius
# ----------------------------------------------------------------------------

def test_grid_neighbor_counts():
    """3x3 grid, pitch 40 m, r_S = 50 m, 10x10 footprints.

    Wall-to-wall distances:
      - cardinal neighbor (pitch 40, w=10): 30 m -> in range
      - diagonal neighbor (pitch sqrt(2)*40, w=10 each side): ~42.4 m -> in range
      - two-apart cardinal: 70 m -> out of range

    So the center has 8 neighbors (4 cardinal + 4 diagonal), corners have 3
    (2 cardinal + 1 diagonal), edge-non-corner buildings have 5 (3 cardinal
    + 2 diagonal).
    """
    gdf = synthetic.grid(n=3, pitch=40.0, width=10.0, length=10.0)
    params = RawMetricParams(r_D=100.0, r_S=50.0, epsilon=0.5, sigma_theta=15.0)
    out = compute_raw_metrics(gdf, params=params, progress=False)

    counts = sorted(out["SS_neighbors"].tolist())
    # 4 corners (3), 4 edge-non-corner (5), 1 center (8)
    assert counts == [3, 3, 3, 3, 5, 5, 5, 5, 8]


# ----------------------------------------------------------------------------
# Orientation sweep — OP monotone in orientation difference (single-pair)
# ----------------------------------------------------------------------------

@pytest.mark.parametrize("angle", [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0])
def test_OP_monotone_decreasing_with_orientation_offset(angle):
    """As the orientation offset grows from 0 to 90 deg, OP_raw should fall
    monotonically (sigma_theta=15) while DP_raw stays at its parallel value
    when spacing is large enough that rotation doesn't change wall distance
    appreciably.

    This is the building-block test for sensitivity sweeps over orientation.
    """
    gdf = synthetic.pair(spacing=50.0, orientation_offset_deg=angle,
                         width=10.0, length=20.0)
    params = RawMetricParams(r_D=100.0, r_S=80.0, epsilon=0.5, sigma_theta=15.0)
    out = compute_raw_metrics(gdf, params=params, progress=False)

    actual_dist = gdf.geometry.iloc[0].distance(gdf.geometry.iloc[1])
    expected_dp = 1.0 / (actual_dist + params.epsilon)
    expected_op = orientation_factor(angle, params.sigma_theta) * expected_dp

    assert out["DP_raw"].iloc[0] == pytest.approx(expected_dp, rel=1e-9)
    assert out["OP_raw"].iloc[0] == pytest.approx(expected_op, abs=1e-12)


# ----------------------------------------------------------------------------
# Centroid coordinates
# ----------------------------------------------------------------------------

def test_centroid_columns_match_rep_points():
    """``cent_x`` / ``cent_y`` should equal the representative-point coords."""
    gdf = synthetic.grid(n=2, pitch=30.0, width=10.0, length=10.0)
    out = compute_raw_metrics(gdf, params=RawMetricParams(), progress=False)
    for i in range(len(out)):
        rep = out.geometry.iloc[i].representative_point()
        assert out["cent_x"].iloc[i] == pytest.approx(rep.x, abs=1e-9)
        assert out["cent_y"].iloc[i] == pytest.approx(rep.y, abs=1e-9)


# ----------------------------------------------------------------------------
# Nearest-neighbor proximity
# ----------------------------------------------------------------------------

def test_nn_proximity_pair_distance_and_bearing():
    """For two rectangles separated by a known wall-to-wall distance, the
    nearest-neighbor distance should equal that spacing, and the bearing
    from B1 (west) to B2 (east) should be 90 deg (compass east); vice versa
    270 deg.
    """
    spacing = 12.0
    gdf = synthetic.pair(spacing=spacing, orientation_offset_deg=0.0,
                         width=10.0, length=20.0)
    out = compute_raw_metrics(gdf, params=RawMetricParams(r_NN=100.0),
                              progress=False)

    assert out["dist_to_nearest_building"].iloc[0] == pytest.approx(spacing, rel=1e-9)
    assert out["dist_to_nearest_building"].iloc[1] == pytest.approx(spacing, rel=1e-9)
    # B2 sits east of B1 -> bearing from B1 to B2 is 90 deg (compass east).
    assert out["bearing_to_nearest_building"].iloc[0] == pytest.approx(90.0, abs=1e-9)
    # B1 sits west of B2 -> bearing from B2 to B1 is 270 deg.
    assert out["bearing_to_nearest_building"].iloc[1] == pytest.approx(270.0, abs=1e-9)


def test_nn_proximity_isolated_building_is_nan():
    """A solo building should have NaN for both NN proximity columns."""
    import math
    gdf = synthetic.isolated_building(width=10.0, length=20.0)
    out = compute_raw_metrics(gdf, params=RawMetricParams(r_NN=200.0),
                              progress=False)
    assert math.isnan(out["dist_to_nearest_building"].iloc[0])
    assert math.isnan(out["bearing_to_nearest_building"].iloc[0])


def test_nn_proximity_out_of_range_is_nan():
    """If the nearest building is beyond ``r_NN``, both columns should be NaN."""
    import math
    spacing = 50.0
    gdf = synthetic.pair(spacing=spacing, orientation_offset_deg=0.0,
                         width=10.0, length=20.0)
    out = compute_raw_metrics(gdf, params=RawMetricParams(r_NN=20.0),
                              progress=False)
    assert math.isnan(out["dist_to_nearest_building"].iloc[0])
    assert math.isnan(out["bearing_to_nearest_building"].iloc[0])
