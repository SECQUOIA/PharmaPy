"""Shortcut distillation design regressions for issue #72."""

import numpy as np
import pytest

from assimulo_helpers import import_module_with_assimulo_stub


pytestmark = pytest.mark.unit


def _import_distillation_module(monkeypatch):
    """Import Distillation while stubbing only optional Assimulo symbols.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest cleanup fixture used by the shared optional-import helper.

    Returns
    -------
    module
        Imported ``PharmaPy.Distillation`` module.
    """
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
    z_feed = np.array([0.40, 0.60])  # [-], feed mole fractions

    x_dist = np.array([0.90, 0.10])  # [-], distillate mole fractions
    x_bottom = np.array([0.10, 0.90])  # [-], bottoms mole fractions
    dist_flowrate = 10.0  # [mol/s]
    bot_flowrate = 15.0  # [mol/s]
    min_reflux = 2.0  # [-], Lmin/D
    num_min = 4.0  # [-], minimum equilibrium-stage count

    class ShortcutColumn(module.DistillationColumn):
        """Shortcut-design double with deterministic correlations."""

        def global_material_bce(self, received_z_feed=None):
            np.testing.assert_allclose(received_z_feed, z_feed)
            return x_dist, x_bottom, dist_flowrate, bot_flowrate

        def calc_underwood_min_reflux(self, *args):
            return min_reflux

        def calc_num_min(self, received_x_dist, received_x_bottom):
            np.testing.assert_allclose(received_x_dist, x_dist)
            np.testing.assert_allclose(received_x_bottom, x_bottom)
            return num_min

    column = ShortcutColumn(
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
    column.z_feed = z_feed

    result = column.calculate_shortcut_design()

    assert result["min_reflux"] == pytest.approx(min_reflux)
    assert result["reflux"] == pytest.approx(3.6)  # 1.8 [-] * Rmin 2.0 [-]


def test_default_shortcut_reflux_ratio_is_pinned(monkeypatch):
    """The documented default multiplier remains 1.5 [-]."""
    module = _import_distillation_module(monkeypatch)
    z_feed = np.array([0.40, 0.60])  # [-], feed mole fractions
    x_dist = np.array([0.90, 0.10])  # [-], distillate mole fractions
    x_bottom = np.array([0.10, 0.90])  # [-], bottoms mole fractions

    class ShortcutColumn(module.DistillationColumn):
        """Shortcut-design double with deterministic minimum reflux."""

        def global_material_bce(self, received_z_feed=None):
            np.testing.assert_allclose(received_z_feed, z_feed)
            dist_flowrate = 10.0  # [mol/s]
            bot_flowrate = 15.0  # [mol/s]
            return x_dist, x_bottom, dist_flowrate, bot_flowrate

        def calc_underwood_min_reflux(self, *args):
            return 2.0  # [-]

        def calc_num_min(self, received_x_dist, received_x_bottom):
            np.testing.assert_allclose(received_x_dist, x_dist)
            np.testing.assert_allclose(received_x_bottom, x_bottom)
            return 4.0  # [-]

    column = ShortcutColumn(
        pres=101325.0,  # [Pa]
        q_feed=1.0,  # [-]
        LK="light",
        HK="heavy",
        perc_LK=95.0,  # [%]
        perc_HK=5.0,  # [%]
        num_plates=8,  # [-]
        num_feed=4,  # [-]
    )
    column.z_feed = z_feed

    result = column.calculate_shortcut_design()

    assert module.DEFAULT_REFLUX_TO_MINIMUM_RATIO == pytest.approx(1.5)
    assert result["reflux"] == pytest.approx(3.0)  # 1.5 [-] * Rmin 2.0 [-]


def test_backward_compatible_shortcut_aliases(monkeypatch):
    module = _import_distillation_module(monkeypatch)
    expected = {
        "material_balances": {},  # flows [mol/s], fractions [-]
        "min_reflux": 1.2,  # [-]
        "num_min": 4.0,  # [-], equilibrium-stage count
        "reflux": 1.8,  # [-]
        "num_plates": 7.0,  # [-], equilibrium-stage count
        "num_feed": 4.0,  # [-], tray count from the top
    }

    class AliasColumn(module.DistillationColumn):
        """Shortcut alias double with deterministic delegate methods."""

        def calculate_shortcut_design(self, time=None):
            assert time is None
            return expected

        def calc_underwood_min_reflux(self, *args):
            return 0.7  # [-]

    column = AliasColumn(
        pres=101325.0,  # [Pa]
        q_feed=1.0,  # [-], saturated liquid feed
        LK="light",
        HK="heavy",
        perc_LK=95.0,  # [%]
        perc_HK=5.0,  # [%]
    )

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
    """Underwood minimum reflux includes feed quality in the root target."""
    module = _import_distillation_module(monkeypatch)
    alpha = np.array([4.0, 1.0])  # [-], relative volatility to HK
    z_feed = np.array([0.40, 0.60])  # [-], feed mole fractions
    x_dist = np.array([0.90, 0.10])  # [-], distillate mole fractions
    x_bottom = np.array([0.05, 0.95])  # [-], bottoms mole fractions
    dist_flowrate = 10.0  # [mol/s]
    bot_flowrate = 15.0  # [mol/s]

    class UnderwoodColumn(module.DistillationColumn):
        """Underwood numeric double with deterministic relative volatility."""

        def get_alpha(self, pres, x_frac):
            """Return the fixture relative-volatility vector.

            Parameters
            ----------
            pres : float
                Column pressure [Pa].
            x_frac : ndarray
                Feed mole fractions [-].

            Returns
            -------
            ndarray
                Relative volatilities referenced to the heavy key [-].
            """
            np.testing.assert_allclose(pres, 101325.0)
            np.testing.assert_allclose(x_frac, z_feed)
            return alpha

    column = UnderwoodColumn(
        pres=101325.0,  # [Pa]
        q_feed=0.8,  # [-], molar liquid fraction in the feed
        LK="light",
        HK="heavy",
        perc_LK=95.0,  # [%]
        perc_HK=5.0,  # [%]
    )
    column.LK_index = 0
    column.HK_index = 1

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


def test_column_startup_accepts_deprecated_time_heuristics_keyword(monkeypatch):
    """The old startup keyword remains available with a deprecation warning."""
    module = _import_distillation_module(monkeypatch)

    class StartupColumn(module.DynamicDistillation):
        """Dynamic column double with deterministic shortcut design."""

        def calculate_shortcut_design(self, time=None):
            self.received_time = time  # [s]
            return {
                "material_balances": {
                    "bottom_flow": 4.0,  # [mol/s]
                    "dist_flow": 6.0,  # [mol/s]
                    "x_dist": np.array([0.9, 0.1]),  # [-]
                    "x_bottom": np.array([0.2, 0.8]),  # [-]
                },
                "min_reflux": 0.7,  # [-]
                "num_min": 3.5,  # [-]
                "reflux": 1.05,  # [-]
                "num_plates": 8.0,  # [-]
                "num_feed": 4.0,  # [-]
            }

    column = StartupColumn(
        pres=101325.0,  # [Pa]
        q_feed=1.0,  # [-]
        LK="light",
        HK="heavy",
        perc_LK=95.0,  # [%]
        perc_HK=5.0,  # [%]
    )

    with pytest.warns(DeprecationWarning, match="time_heuristics is deprecated"):
        column.column_startup(time_heuristics=12.0)

    assert column.received_time == pytest.approx(12.0)
    assert column.num_plates == 8
    assert column.num_feed == 4
    assert column.reflux == pytest.approx(1.05)
