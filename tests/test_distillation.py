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
        def __init__(self, design):
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
            return self.design

        def calc_plates(self, **kwargs):
            self.calc_plates_calls.append(kwargs)

        def retrieve_results(self, *args, **kwargs):
            raise AssertionError("solve_unit should not retrieve results twice")

    column = ShortcutDesignColumn(expected)

    result = column.solve_unit(solve_ss=True)

    assert result == expected
    assert column.calc_plates_calls == [{
        "x_dist": [0.95, 0.05],
        "x_bottom": [0.05, 0.95],
        "dist_flow": 1.0,
        "bottom_flow": 2.0,
        "reflux": 1.8,
        "num_plates": 7,
    }]
