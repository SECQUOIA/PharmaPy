"""Regression tests for DynamicExtractor default equilibrium routing."""

import importlib
import sys
import types

import numpy as np
import pytest


pytestmark = pytest.mark.unit


def _load_dynamic_extraction(monkeypatch):
    """Import DynamicExtraction with a minimal solver stub for unit tests."""
    assimulo = types.ModuleType("assimulo")
    problem = types.ModuleType("assimulo.problem")
    solvers = types.ModuleType("assimulo.solvers")

    class _UnavailableSolver:
        def __init__(self, *args, **kwargs):
            raise RuntimeError("solver stub is not executable")

    problem.Implicit_Problem = object
    solvers.IDA = _UnavailableSolver
    solvers.Radau5DAE = _UnavailableSolver

    monkeypatch.setitem(sys.modules, "assimulo", assimulo)
    monkeypatch.setitem(sys.modules, "assimulo.problem", problem)
    monkeypatch.setitem(sys.modules, "assimulo.solvers", solvers)
    monkeypatch.delitem(sys.modules, "PharmaPy.DynamicExtraction", raising=False)

    return importlib.import_module("PharmaPy.DynamicExtraction")


class _ActivityPhase:
    def __init__(self):
        self.calls = []

    def getActivityCoeff(self, method, mole_frac, temp):
        self.calls.append((method, np.asarray(mole_frac), np.asarray(temp)))
        return 1.0 + 0.5 * np.asarray(mole_frac)  # gamma_i [-]


class _StopAfterConstruction(Exception):
    """Raised by the stub to stop initialization at a known point."""


def test_default_k_fun_uses_selected_activity_model(monkeypatch):
    """The default DynamicExtractor k_fun returns gamma_light/gamma_heavy."""
    module = _load_dynamic_extraction(monkeypatch)
    extractor = module.DynamicExtractor(num_stages=2, gamma_model="UNIFAC")
    extractor.Liquid_1 = _ActivityPhase()

    x_light = np.array([[0.2, 0.8], [0.4, 0.6]])  # [-]
    x_heavy = np.array([[0.5, 0.5], [0.1, 0.9]])  # [-]
    temp = np.array([298.15, 301.15])  # K

    k_i = extractor.k_fun(x_light, x_heavy, temp)  # K_i [-]

    expected = (1.0 + 0.5 * x_light) / (1.0 + 0.5 * x_heavy)  # [-]
    assert k_i == pytest.approx(expected)
    assert [call[0] for call in extractor.Liquid_1.calls] == [
        "UNIFAC",
        "UNIFAC",
    ]
    assert extractor.Liquid_1.calls[0][2] == pytest.approx(temp)


def test_default_k_fun_rejects_unknown_activity_model(monkeypatch):
    """DynamicExtractor rejects unknown default activity-model selectors."""
    module = _load_dynamic_extraction(monkeypatch)

    with pytest.raises(ValueError, match="gamma_model must be one of"):
        module.DynamicExtractor(num_stages=1, gamma_model="uniquac")


def test_initialize_model_passes_default_k_fun_to_batch_extractor(monkeypatch):
    """The default constructor supplies a callable equilibrium function."""
    module = _load_dynamic_extraction(monkeypatch)
    captured = {}

    class _BatchExtractor:
        def __init__(self, k_fun=None):
            captured["k_fun"] = k_fun

        def solve_unit(self):
            raise _StopAfterConstruction

    monkeypatch.setattr(module, "BatchExtractor", _BatchExtractor)

    extractor = module.DynamicExtractor(num_stages=1)
    extractor.Liquid_1 = object()

    with pytest.raises(_StopAfterConstruction):
        extractor.initialize_model()

    assert captured["k_fun"] is extractor.k_fun
    assert callable(captured["k_fun"])
