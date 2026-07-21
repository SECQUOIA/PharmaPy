"""Regression tests for DynamicExtractor default equilibrium routing."""

import numpy as np
import pytest

from assimulo_helpers import import_module_with_assimulo_stub


pytestmark = pytest.mark.unit


def _load_dynamic_extraction(monkeypatch):
    """Import DynamicExtraction while stubbing only optional Assimulo symbols.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest cleanup fixture used by the shared optional-import helper.

    Returns
    -------
    module
        Imported ``PharmaPy.DynamicExtraction`` module.
    """
    return import_module_with_assimulo_stub(
        monkeypatch,
        "PharmaPy.DynamicExtraction",
        solvers={"IDA": object, "Radau5DAE": object},
        problem={"Implicit_Problem": object},
    )


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
    temp = np.array([298.15, 301.15])  # [K]

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

    # Lowercase ``uniquac`` is deliberate: activity-model selectors are
    # case-sensitive and must match the exact Phases.getActivityCoeff branch.
    with pytest.raises(ValueError, match="gamma_model must be one of"):
        module.DynamicExtractor(num_stages=1, gamma_model="uniquac")


def test_initialize_model_passes_default_k_fun_to_batch_extractor(monkeypatch):
    """The default constructor supplies ``BatchExtractor`` the selected model.

    This pins the DynamicExtractor -> BatchExtractor handoff without running a
    full static extraction solve. The real solve depends on thermodynamic phase
    objects, outlet construction, and broader stage-wise K_i behavior tracked
    separately in #123; the sentinel stops at the fixed boundary.
    """
    module = _load_dynamic_extraction(monkeypatch)
    captured = {}

    class _BatchExtractor:
        def __init__(self, k_fun=None, gamma_method="UNIQUAC"):
            captured["k_fun"] = k_fun
            captured["gamma_method"] = gamma_method

        def solve_unit(self):
            raise _StopAfterConstruction

    monkeypatch.setattr(module, "BatchExtractor", _BatchExtractor)

    extractor = module.DynamicExtractor(num_stages=1, gamma_model="ideal")
    extractor.Liquid_1 = object()

    with pytest.raises(_StopAfterConstruction):
        extractor.initialize_model()

    assert captured["k_fun"] is extractor.k_fun
    assert callable(captured["k_fun"])
    assert captured["gamma_method"] == "ideal"
