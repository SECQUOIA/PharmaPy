"""Shortcut distillation design regressions for issue #72."""

import numpy as np
import pytest

from assimulo_helpers import import_module_with_assimulo_stub


pytestmark = pytest.mark.unit


def _import_distillation_module(monkeypatch):
    module = import_module_with_assimulo_stub(
        monkeypatch,
        "PharmaPy.Distillation",
        solvers={"IDA": object},
        problem={"Implicit_Problem": object},
    )
    return module


def test_positive_reflux_below_minimum_uses_configured_shortcut_ratio(
        monkeypatch):
    module = _import_distillation_module(monkeypatch)
    column = module.DistillationColumn(
        pres=101325.0,  # [Pa]
        q_feed=1.0,  # [-], saturated liquid feed
        LK="light",
        HK="heavy",
        perc_LK=95.0,  # [%]
        perc_HK=5.0,  # [%]
        reflux=1.0,  # [-], below the computed minimum
        num_plates=8,  # [-], equilibrium-stage count
        num_feed=4,  # [-], tray count from the top
        reflux_to_minimum_ratio=1.8,  # [-], operating reflux/Rmin
    )
    z_feed = np.array([0.40, 0.60])  # [-], feed mole fractions
    column.z_feed = z_feed

    x_dist = np.array([0.90, 0.10])  # [-], distillate mole fractions
    x_bottom = np.array([0.10, 0.90])  # [-], bottoms mole fractions
    dist_flowrate = 10.0  # [mol/s]
    bot_flowrate = 15.0  # [mol/s]
    min_reflux = 2.0  # [-], Lmin/D
    num_min = 4.0  # [-], minimum equilibrium-stage count

    def fake_global_material_bce(received_z_feed):
        np.testing.assert_allclose(received_z_feed, z_feed)
        return x_dist, x_bottom, dist_flowrate, bot_flowrate

    def fake_calc_min_reflux(*args):
        return min_reflux

    def fake_calc_num_min(received_x_dist, received_x_bottom):
        np.testing.assert_allclose(received_x_dist, x_dist)
        np.testing.assert_allclose(received_x_bottom, x_bottom)
        return num_min

    monkeypatch.setattr(column, "global_material_bce", fake_global_material_bce)
    monkeypatch.setattr(
        column, "calc_underwood_min_reflux", fake_calc_min_reflux)
    monkeypatch.setattr(column, "calc_num_min", fake_calc_num_min)

    result = column.calculate_shortcut_design()

    assert result["min_reflux"] == pytest.approx(min_reflux)
    assert result["reflux"] == pytest.approx(
        column.reflux_to_minimum_ratio * min_reflux)


def test_backward_compatible_shortcut_aliases(monkeypatch):
    module = _import_distillation_module(monkeypatch)
    column = module.DistillationColumn(
        pres=101325.0,  # [Pa]
        q_feed=1.0,  # [-], saturated liquid feed
        LK="light",
        HK="heavy",
        perc_LK=95.0,  # [%]
        perc_HK=5.0,  # [%]
    )
    expected = {
        "material_balances": {},  # [mol/s, -]
        "min_reflux": 1.2,  # [-]
        "num_min": 4.0,  # [-], equilibrium-stage count
        "reflux": 1.8,  # [-]
        "num_plates": 7.0,  # [-], equilibrium-stage count
        "num_feed": 4.0,  # [-], tray count from the top
    }

    def fake_calculate_shortcut_design(time=None):
        assert time is None
        return expected

    def fake_calc_underwood_min_reflux(*args):
        return 0.7  # [-]

    monkeypatch.setattr(
        column, "calculate_shortcut_design", fake_calculate_shortcut_design)
    monkeypatch.setattr(
        column, "calc_underwood_min_reflux", fake_calc_underwood_min_reflux)

    with pytest.warns(
            DeprecationWarning,
            match="calculate_heuristics is deprecated"):
        assert column.calculate_heuristics() is expected

    with pytest.warns(
            DeprecationWarning,
            match="calc_min_reflux is deprecated"):
        assert column.calc_min_reflux(None, None, None, None) == pytest.approx(
            0.7)


def test_underwood_min_reflux_includes_feed_quality_in_target(monkeypatch):
    module = _import_distillation_module(monkeypatch)
    column = module.DistillationColumn(
        pres=101325.0,  # [Pa]
        q_feed=0.8,  # [-], molar liquid fraction in the feed
        LK="light",
        HK="heavy",
        perc_LK=95.0,  # [%]
        perc_HK=5.0,  # [%]
    )
    column.LK_index = 0
    column.HK_index = 1

    alpha = np.array([4.0, 1.0])  # [-], relative volatility to HK
    z_feed = np.array([0.40, 0.60])  # [-], feed mole fractions
    x_dist = np.array([0.90, 0.10])  # [-], distillate mole fractions
    x_bottom = np.array([0.05, 0.95])  # [-], bottoms mole fractions
    dist_flowrate = 10.0  # [mol/s]
    bot_flowrate = 15.0  # [mol/s]

    def fake_get_alpha(pres, x_frac):
        assert pres == pytest.approx(101325.0)
        np.testing.assert_allclose(x_frac, z_feed)
        return alpha

    monkeypatch.setattr(column, "get_alpha", fake_get_alpha)

    min_reflux = column.calc_underwood_min_reflux(
        x_dist=x_dist,
        x_bot=x_bottom,
        dist_flowrate=dist_flowrate,
        bot_flowrate=bot_flowrate,
        z_feed=z_feed,
    )

    # With q_feed = 0.8 [-], the first Underwood target is 1-q = 0.2 [-].
    # For this fixture phi = 2.0 [-], Vmin = 17.0 [mol/s], and Rmin = 0.7 [-].
    assert min_reflux == pytest.approx(0.7)
