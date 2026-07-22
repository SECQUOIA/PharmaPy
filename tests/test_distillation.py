"""Distillation solve-unit regressions."""

import pytest

pytest.importorskip("assimulo")
from PharmaPy.Distillation import DistillationColumn


pytestmark = [pytest.mark.assimulo, pytest.mark.unit]


def test_steady_state_solve_uses_calc_plates_results():
    """``solve_unit`` forwards shortcut-design results to ``calc_plates``.

    A small subclass test double avoids monkeypatching the instance while
    keeping the assertion focused on the ``solve_unit`` orchestration branch.
    The shortcut-design and plate-calculation numerics are covered by their
    dedicated tests, so this test stops before result retrieval.
    """
    material_balances = {
        "x_dist": [0.95, 0.05],  # [-]
        "x_bottom": [0.05, 0.95],  # [-]
        "dist_flow": 1.0,  # [mol/s]
        "bottom_flow": 2.0,  # [mol/s]
    }
    expected = {
        "material_balances": material_balances,
        "min_reflux": 1.2,  # [-]
        "num_min": 4,  # [-]
        "reflux": 1.8,  # [-]
        "num_plates": 7,  # [-]
    }

    class ShortcutDesignColumn(DistillationColumn):
        """Column test double for shortcut-design handoff."""

        def __init__(self, design):
            """Store a deterministic shortcut design.

            Parameters
            ----------
            design : dict
                Shortcut-design result. Mole fractions are [-], molar flows are
                [mol/s], reflux values are [-], and stage counts are [-].
            """
            super().__init__(
                pres=101325.0,  # [Pa]
                q_feed=1.0,  # [-]
                LK="ethanol",
                HK="water",
                perc_LK=95.0,  # [%]
                perc_HK=5.0,  # [%]
            )
            self.design = design
            self.calc_plates_calls = []

        def calculate_shortcut_design(self):
            """Return the deterministic shortcut-design result.

            Returns
            -------
            dict
                Shortcut-design result with mole fractions [-], molar flows
                [mol/s], reflux values [-], and stage counts [-].
            """
            return self.design

        def calc_plates(self, **kwargs):
            """Capture the plate-calculation keyword arguments.

            Parameters
            ----------
            **kwargs
                Plate-calculation inputs. Mole fractions are [-], molar flows
                are [mol/s], reflux is [-], and stage counts are [-].

            Returns
            -------
            None
                The call is recorded in ``calc_plates_calls``.
            """
            self.calc_plates_calls.append(kwargs)

        def retrieve_results(self, *args, **kwargs):
            """Prevent accidental result retrieval in this handoff test.

            Parameters
            ----------
            *args
                Positional result-retrieval arguments.
            **kwargs
                Keyword result-retrieval arguments.

            Raises
            ------
            AssertionError
                Always, because this test stops before retrieval.
            """
            raise AssertionError("solve_unit should not retrieve results twice")

    column = ShortcutDesignColumn(expected)

    result = column.solve_unit(solve_ss=True)

    assert result == expected
    assert column.calc_plates_calls == [{
        "x_dist": [0.95, 0.05],  # [-]
        "x_bottom": [0.05, 0.95],  # [-]
        "dist_flow": 1.0,  # [mol/s]
        "bottom_flow": 2.0,  # [mol/s]
        "reflux": 1.8,  # [-]
        "num_plates": 7,  # [-]
    }]
