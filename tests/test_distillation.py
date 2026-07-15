import pytest

pytest.importorskip("assimulo")
from PharmaPy.Distillation import DistillationColumn


pytestmark = [pytest.mark.assimulo, pytest.mark.unit]


def test_steady_state_solve_uses_calc_plates_results(monkeypatch):
    column = DistillationColumn(
        pres=101325,
        q_feed=1.0,
        LK="ethanol",
        HK="water",
        perc_LK=95.0,
        perc_HK=5.0,
    )
    material_balances = {
        "x_dist": [0.95, 0.05],
        "x_bottom": [0.05, 0.95],
        "dist_flow": 1.0,
        "bottom_flow": 2.0,
    }
    expected = {
        "material_balances": material_balances,
        "min_reflux": 1.2,
        "num_min": 4,
        "reflux": 1.8,
        "num_plates": 7,
    }
    calls = []

    def fake_calculate_heuristics():
        return expected

    def fake_calc_plates(**kwargs):
        calls.append(kwargs)

    def fail_retrieve_results(*args, **kwargs):
        raise AssertionError("solve_unit should not retrieve results twice")

    monkeypatch.setattr(column, "calculate_heuristics", fake_calculate_heuristics)
    monkeypatch.setattr(column, "calc_plates", fake_calc_plates)
    monkeypatch.setattr(column, "retrieve_results", fail_retrieve_results)

    result = column.solve_unit(solve_ss=True)

    assert result == expected
    assert calls == [{
        "x_dist": [0.95, 0.05],
        "x_bottom": [0.05, 0.95],
        "dist_flow": 1.0,
        "bottom_flow": 2.0,
        "reflux": 1.8,
        "num_plates": 7,
    }]
