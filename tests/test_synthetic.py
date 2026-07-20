"""Accuracy tests against analytically known synthetic geometries.

Each test builds a small fixture from :mod:`ssdd.synthetic`, runs
:func:`ssdd.pipeline.compute_raw_metrics`, and checks ``SD`` / ``SS``
against hand-computed expected values.

The math behind each expectation is spelled out in the docstring of the test
so it doubles as documentation of how the metrics behave on simple geometry.
"""

from __future__ import annotations

import math

import pytest

from ssdd import synthetic
from ssdd.geometry import angle_difference_deg, kernel_value
from ssdd.pipeline import RawMetricParams, compute_raw_metrics


# ----------------------------------------------------------------------------
# Kernel family
# ----------------------------------------------------------------------------

def test_kernel_family_values_at_half():
    """K_n(0.5) = (1 - 0.25)^n and truncation outside [0, 1]."""
    assert kernel_value(0.5, "uniform") == 1.0
    assert kernel_value(0.5, "epanechnikov") == pytest.approx(0.75)
    assert kernel_value(0.5, "quartic") == pytest.approx(0.5625)
    assert kernel_value(0.5, "triweight") == pytest.approx(0.421875)
    assert kernel_value(1.0001, "uniform") == 0.0
    with pytest.raises(ValueError):
        kernel_value(0.5, "gaussian")


# ----------------------------------------------------------------------------
# Isolated building
# ----------------------------------------------------------------------------

def test_isolated_building_scores_zero():
    """One 10x20 rectangle, no neighbors: SD = 0 and SS = 0 (isolated)."""
    gdf = synthetic.isolated_building(width=10.0, length=20.0)
    out = compute_raw_metrics(gdf, params=RawMetricParams(), progress=False)

    assert out["SD"].iloc[0] == 0.0
    assert out["SS"].iloc[0] == 0.0
    assert out["bld_area"].iloc[0] == pytest.approx(200.0)


# ----------------------------------------------------------------------------
# Parallel pair — both metrics recover closed forms
# ----------------------------------------------------------------------------

def test_pair_recovers_closed_forms():
    """Two 10x20 rectangles (area 200), walls 10 m apart, axis-aligned.

    Expected with defaults (uniform kernel, root_area weight, r_D=200, eps=0.5):
      SS = 1 / (10 + 0.5)
      centroid distance = spacing + width = 20 m <= r_D, K_uniform = 1
      SD = sqrt(200) / (pi * 200^2)
    """
    spacing, width, length = 10.0, 10.0, 20.0
    gdf = synthetic.pair(spacing=spacing, orientation_offset_deg=0.0,
                         width=width, length=length)
    p = RawMetricParams()
    out = compute_raw_metrics(gdf, params=p, progress=False)

    expected_ss = 1.0 / (spacing + p.epsilon)
    assert out["SS"].iloc[0] == pytest.approx(expected_ss, rel=1e-9)
    assert out["SS"].iloc[1] == pytest.approx(expected_ss, rel=1e-9)

    expected_sd = math.sqrt(width * length) / (math.pi * p.r_D ** 2)
    assert out["SD"].iloc[0] == pytest.approx(expected_sd, rel=1e-9)
    assert out["SD"].iloc[1] == pytest.approx(expected_sd, rel=1e-9)


def test_pair_weight_variants():
    """Same pair; SD scales as w_j: unit -> 1, area -> 200, root_area -> sqrt(200)."""
    gdf = synthetic.pair(spacing=10.0, orientation_offset_deg=0.0,
                         width=10.0, length=20.0)
    base = math.pi * 200.0 ** 2
    for weight, w in [("unit", 1.0), ("area", 200.0), ("root_area", math.sqrt(200.0))]:
        out = compute_raw_metrics(
            gdf, params=RawMetricParams(weight=weight), progress=False)
        assert out["SD"].iloc[0] == pytest.approx(w / base, rel=1e-9), weight


def test_pair_kernel_variants():
    """Pair at centroid distance 20 m with r_D=100 -> u=0.2; SD scales as K(u)."""
    gdf = synthetic.pair(spacing=10.0, orientation_offset_deg=0.0,
                         width=10.0, length=20.0)
    u = 0.2
    for kernel in ("uniform", "epanechnikov", "quartic", "triweight"):
        out = compute_raw_metrics(
            gdf, params=RawMetricParams(r_D=100.0, kernel=kernel, weight="unit"),
            progress=False)
        expected = kernel_value(u, kernel) / (math.pi * 100.0 ** 2)
        assert out["SD"].iloc[0] == pytest.approx(expected, rel=1e-9), kernel


# ----------------------------------------------------------------------------
# SS truncation and saturation
# ----------------------------------------------------------------------------

def test_touching_pair_saturates_at_one_over_epsilon():
    """Two rectangles sharing a wall: d = 0, so SS = 1 / epsilon."""
    eps = 0.5
    gdf = synthetic.touching_pair(orientation_offset_deg=0.0)
    out = compute_raw_metrics(gdf, params=RawMetricParams(epsilon=eps),
                              progress=False)
    assert out["SS"].iloc[0] == pytest.approx(1.0 / eps, rel=1e-9)
    assert out["SS"].iloc[1] == pytest.approx(1.0 / eps, rel=1e-9)


def test_ss_nn_ignores_r_s():
    """nn uses the TRUE nearest neighbour and ignores r_S: a pair 50 m apart
    still scores 1/(50+eps) even when r_S=20 (there is no isolation cutoff)."""
    gdf = synthetic.pair(spacing=50.0, orientation_offset_deg=0.0,
                         width=10.0, length=20.0)
    p = RawMetricParams(r_S=20.0)
    out = compute_raw_metrics(gdf, params=p, progress=False)
    d = gdf.geometry.iloc[0].distance(gdf.geometry.iloc[1])
    assert out["SS"].iloc[0] == pytest.approx(1.0 / (d + p.epsilon), rel=1e-9)


def test_ss_averaging_still_respects_r_s():
    """The averaging aggs DO truncate at r_S: with r_S=20 the 25 m neighbour is
    dropped, leaving only the 10 m one, so power1 = 1/(10+eps)."""
    specs = [(0.0, 0.0, 10.0, 10.0),
             (20.0, 0.0, 10.0, 10.0),   # wall 10  (kept)
             (0.0, 35.0, 10.0, 10.0)]   # wall 25  (dropped at r_S=20)
    gdf = synthetic.sized_cluster(specs)
    eps = 0.5
    out = compute_raw_metrics(
        gdf, params=RawMetricParams(r_S=20.0, epsilon=eps, agg="power1"),
        progress=False)
    assert out["SS"].iloc[0] == pytest.approx(1.0 / (10.0 + eps), rel=1e-9)


def test_ss_orientation_invariance():
    """SS ignores orientation: rotating the neighbour changes SS only through
    the (actual) wall-to-wall distance, never through an angle weight."""
    p = RawMetricParams()
    for angle in (0.0, 45.0, 90.0):
        gdf = synthetic.pair(spacing=30.0, orientation_offset_deg=angle,
                             width=10.0, length=20.0)
        out = compute_raw_metrics(gdf, params=p, progress=False)
        actual = gdf.geometry.iloc[0].distance(gdf.geometry.iloc[1])
        assert out["SS"].iloc[0] == pytest.approx(
            1.0 / (actual + p.epsilon), rel=1e-9), angle


# ----------------------------------------------------------------------------
# 3x3 grid — SD truncation counts, SS uniform nearest distance
# ----------------------------------------------------------------------------

def test_grid_sd_truncation_counts():
    """3x3 grid, pitch 40 m, 10x10 footprints, r_D = 50, unit weight, uniform K.

    Centroid distances: cardinal 40 (in), diagonal ~56.6 (out), two-apart 80
    (out). So SD * pi * r_D^2 = number of cardinal neighbors: corners 2,
    edges 3, center 4.
    """
    gdf = synthetic.grid(n=3, pitch=40.0, width=10.0, length=10.0)
    out = compute_raw_metrics(
        gdf, params=RawMetricParams(r_D=50.0, weight="unit"), progress=False)
    counts = sorted(round(v * math.pi * 50.0 ** 2) for v in out["SD"])
    assert counts == [2, 2, 2, 2, 3, 3, 3, 3, 4]


def test_grid_ss_uniform_nearest():
    """Same grid: every building's nearest wall is a cardinal neighbor at
    pitch - width = 30 m, so SS = 1 / 30.5 everywhere."""
    gdf = synthetic.grid(n=3, pitch=40.0, width=10.0, length=10.0)
    p = RawMetricParams()
    out = compute_raw_metrics(gdf, params=p, progress=False)
    expected = 1.0 / (30.0 + p.epsilon)
    assert out["SS"].tolist() == pytest.approx([expected] * 9, rel=1e-9)


# ----------------------------------------------------------------------------
# Attributes
# ----------------------------------------------------------------------------

def test_centroid_columns_match_true_centroids():
    """``cent_x`` / ``cent_y`` are true polygon centroids (also SD's rep points)."""
    gdf = synthetic.grid(n=2, pitch=30.0, width=10.0, length=10.0)
    out = compute_raw_metrics(gdf, params=RawMetricParams(), progress=False)
    for i in range(len(out)):
        c = out.geometry.iloc[i].centroid
        assert out["cent_x"].iloc[i] == pytest.approx(c.x, abs=1e-9)
        assert out["cent_y"].iloc[i] == pytest.approx(c.y, abs=1e-9)


def test_phi_deg_attribute_and_angle_folding():
    """phi_deg tracks the long axis; angle_difference_deg folds to [0, 90]."""
    gdf = synthetic.pair(spacing=30.0, orientation_offset_deg=90.0,
                         width=10.0, length=20.0)
    out = compute_raw_metrics(gdf, params=RawMetricParams(), progress=False)
    d = angle_difference_deg(out["phi_deg"].iloc[0], out["phi_deg"].iloc[1])
    assert d == pytest.approx(90.0, abs=1e-6)
    assert angle_difference_deg(10.0, 170.0) == pytest.approx(20.0)


# ----------------------------------------------------------------------------
# Heterogeneous sizes — the SD weighting choices must diverge
# ----------------------------------------------------------------------------

def test_sized_cluster_weight_scaling():
    """Focal at origin with two neighbours of DIFFERENT area, both at centroid
    distance 40 <= r_D=100 (uniform K=1). SD * (pi r_D^2) is the summed weight:
      unit      -> 1 + 1     = 2
      area      -> 100 + 400 = 500
      root_area -> 10 + 20   = 30
    """
    specs = [(0.0, 0.0, 10.0, 10.0),    # focal, area 100
             (40.0, 0.0, 10.0, 10.0),   # neighbour, area 100, dist 40
             (0.0, 40.0, 20.0, 20.0)]   # neighbour, area 400, dist 40
    gdf = synthetic.sized_cluster(specs)
    norm = math.pi * 100.0 ** 2
    for weight, expected in [("unit", 2.0), ("area", 500.0), ("root_area", 30.0)]:
        out = compute_raw_metrics(
            gdf, params=RawMetricParams(r_D=100.0, weight=weight), progress=False)
        assert out["SD"].iloc[0] == pytest.approx(expected / norm, rel=1e-9), weight


# ----------------------------------------------------------------------------
# SS aggregations and orientation weights (the non-default forms)
# ----------------------------------------------------------------------------

def test_ss_agg_power1_averages_inverse_distance():
    """Focal with neighbours at wall distances 10 and 25; agg='power1' averages
    1/(d+eps) over both, while the default nn keeps only the nearest."""
    specs = [(0.0, 0.0, 10.0, 10.0),
             (20.0, 0.0, 10.0, 10.0),   # wall 10
             (0.0, 35.0, 10.0, 10.0)]   # wall 25
    gdf = synthetic.sized_cluster(specs)
    eps = 0.5
    out = compute_raw_metrics(
        gdf, params=RawMetricParams(r_S=50.0, epsilon=eps, agg="power1"),
        progress=False)
    expected = (1.0 / (10.0 + eps) + 1.0 / (25.0 + eps)) / 2.0
    assert out["SS"].iloc[0] == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize("orient", ["cos2", "cos4", "gauss"])
def test_ss_orientation_weight_applied_at_nn(orient):
    """With an orientation kernel, nn SS = g(theta*) / (d* + eps), where theta*
    is the folded angle difference to the nearest neighbour."""
    gdf = synthetic.pair(spacing=30.0, orientation_offset_deg=40.0,
                         width=10.0, length=20.0)
    sigma, eps = 12.0, 0.5
    out = compute_raw_metrics(
        gdf, params=RawMetricParams(r_S=60.0, epsilon=eps, orient=orient, sigma=sigma),
        progress=False)
    theta = angle_difference_deg(out["phi_deg"].iloc[0], out["phi_deg"].iloc[1])
    d = gdf.geometry.iloc[0].distance(gdf.geometry.iloc[1])
    g = {"cos2": math.cos(math.radians(theta)) ** 2,
         "cos4": math.cos(math.radians(theta)) ** 4,
         "gauss": math.exp(-((theta / sigma) ** 2))}[orient]
    assert out["SS"].iloc[0] == pytest.approx(g / (d + eps), rel=1e-9)


# ----------------------------------------------------------------------------
# Irregular (L-shaped) footprint
# ----------------------------------------------------------------------------

def test_lshape_area_and_orientation_equivariance():
    """L area = arm_width*(2*arm_length - arm_width); rotating the footprint by
    delta rotates phi_deg by delta (mod 180)."""
    al, aw = 30.0, 12.0
    out0 = compute_raw_metrics(synthetic.lshape(al, aw, angle_deg=0.0),
                               params=RawMetricParams(), progress=False)
    assert out0["bld_area"].iloc[0] == pytest.approx(aw * (2 * al - aw), rel=1e-9)

    delta = 35.0
    out1 = compute_raw_metrics(synthetic.lshape(al, aw, angle_deg=delta),
                               params=RawMetricParams(), progress=False)
    phi0, phi1 = out0["phi_deg"].iloc[0], out1["phi_deg"].iloc[0]
    assert angle_difference_deg(phi1, (phi0 + delta) % 180.0) == pytest.approx(0.0, abs=1e-6)
