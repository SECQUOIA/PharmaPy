import json
import sys
from types import ModuleType

import numpy as np
import pytest

from PharmaPy.StatsModule import StatisticsClass
from PharmaPy.ThermoModule import ParseDatabase


pytestmark = pytest.mark.unit


def _stub_assimulo_modules(monkeypatch):
    assimulo = ModuleType("assimulo")

    solvers = ModuleType("assimulo.solvers")
    solvers.CVode = object

    problem = ModuleType("assimulo.problem")
    problem.Explicit_Problem = object

    exception = ModuleType("assimulo.exception")
    exception.TerminateSimulation = Exception

    monkeypatch.setitem(sys.modules, "assimulo", assimulo)
    monkeypatch.setitem(sys.modules, "assimulo.solvers", solvers)
    monkeypatch.setitem(sys.modules, "assimulo.problem", problem)
    monkeypatch.setitem(sys.modules, "assimulo.exception", exception)


def _import_batch_cryst(monkeypatch):
    try:
        from PharmaPy.Crystallizers import BatchCryst
    except ModuleNotFoundError as exc:
        if exc.name != "assimulo":
            raise
        _stub_assimulo_modules(monkeypatch)

        from PharmaPy.Crystallizers import BatchCryst

    return BatchCryst


class _FakeEstimationInstance:
    def __init__(self, optimize_fn):
        self.num_params = 2
        self.opt_method = "LM"
        self.optim_options = {}
        self.y_data = None
        self.optimize_fn = optimize_fn


def _make_statistics(optimize_fn):
    stats = StatisticsClass.__new__(StatisticsClass)
    stats.inst = _FakeEstimationInstance(optimize_fn)
    stats.get_bootsamples = lambda num_samples: [np.zeros(num_samples)]

    return stats


def test_parse_database_converts_numeric_fields_to_float_arrays(tmp_path):
    database = {
        "water": {"mw": 18.02, "cas": "7732-18-5"},
        "ethanol": {"mw": 46.07, "cas": "64-17-5"},
    }
    path = tmp_path / "compounds.json"
    path.write_text(json.dumps(database))

    parsed = ParseDatabase(str(path))

    assert isinstance(parsed["mw"], np.ndarray)
    assert parsed["mw"].dtype == np.float64
    np.testing.assert_allclose(sorted(parsed["mw"]), [18.02, 46.07])


def test_parse_database_keeps_non_numeric_fields_as_lists(tmp_path):
    database = {
        "water": {"mw": 18.02, "cas": "7732-18-5"},
        "ethanol": {"mw": 46.07, "cas": "64-17-5"},
    }
    path = tmp_path / "compounds.json"
    path.write_text(json.dumps(database))

    parsed = ParseDatabase(str(path))

    assert isinstance(parsed["cas"], list)
    assert sorted(parsed["cas"]) == ["64-17-5", "7732-18-5"]


def test_bootstrap_params_records_nan_rows_for_failed_optimizations():
    def failing_optimize(**kwargs):
        raise RuntimeError("solver blew up")

    stats = _make_statistics(failing_optimize)

    stats.bootstrap_params(num_samples=3)

    assert stats.boot_params.shape == (3, 2)
    assert np.isnan(stats.boot_params).all()


def test_bootstrap_params_does_not_swallow_keyboard_interrupt():
    def interrupted_optimize(**kwargs):
        raise KeyboardInterrupt

    stats = _make_statistics(interrupted_optimize)

    with pytest.raises(KeyboardInterrupt):
        stats.bootstrap_params(num_samples=3)


def test_batch_cryst_warns_when_ad_requested_without_jax(monkeypatch):
    BatchCryst = _import_batch_cryst(monkeypatch)

    monkeypatch.setitem(sys.modules, "jax", None)
    monkeypatch.delitem(sys.modules, "jax.numpy", raising=False)

    with pytest.warns(RuntimeWarning, match="jax"):
        crystallizer = BatchCryst(target_comp="solute", jac_type="AD")

    assert crystallizer.jac_type == "AD"
    assert type(crystallizer).np is np
