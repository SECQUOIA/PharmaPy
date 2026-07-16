"""Regression tests for extractor operating-mode flow/amount routing.

Continuous extraction should consume and return mole-flow rates [mol/s].
Batch extraction should keep material amounts [mol].
"""

import numpy as np
import pytest

from PharmaPy.Extractors import BatchExtractor, ContinuousExtractor


pytestmark = pytest.mark.unit


class _Inlet:
    name_species = ["solute", "solvent"]
    mole_flow = 8.0  # mol/s
    moles = 3.0  # mol
    temp = 298.15  # K
    pres = 101325.0  # Pa
    mole_frac = np.array([0.25, 0.75])  # [-]
    path_data = "dummy-path"


class _DummyOutlet:
    def __init__(self, path, mole_frac, temp, pres, **amount):
        self.path = path
        self.mole_frac = mole_frac  # [-]
        self.temp = temp  # K
        self.pres = pres  # Pa
        self.amount = amount  # mole_flow [mol/s] or moles [mol]

    def getDensity(self, basis):
        assert basis == "mole"
        return 1000.0  # arbitrary molar density [mol/m^3]


class _DummyStream(_DummyOutlet):
    pass


class _DummyPhase(_DummyOutlet):
    pass


def test_continuous_extractor_uses_mole_flow_and_stream_outlets(monkeypatch):
    """Continuous extractor routes inlet and outlet quantities as mol/s."""
    created = {"stream": [], "phase": []}

    def stream_factory(*args, **kwargs):
        outlet = _DummyStream(*args, **kwargs)
        created["stream"].append(outlet)
        return outlet

    def phase_factory(*args, **kwargs):
        outlet = _DummyPhase(*args, **kwargs)
        created["phase"].append(outlet)
        return outlet

    monkeypatch.setattr("PharmaPy.Extractors.LiquidStream", stream_factory)
    monkeypatch.setattr("PharmaPy.Extractors.LiquidPhase", phase_factory)

    extractor = ContinuousExtractor()
    extractor.Inlet = _Inlet()

    assert extractor.oper_mode == "Continuous"
    assert extractor.in_flow == pytest.approx(_Inlet.mole_flow)

    extractor.retrieve_results((
        0.25,  # phase fraction [-]
        np.array([0.4, 0.6]),  # liquid-a mole fractions [-]
        np.array([0.1, 0.9]),  # liquid-b mole fractions [-]
        {"error": 0.0, "num_iter": 1},
    ))

    assert not created["phase"]
    assert len(created["stream"]) == 2
    assert {tuple(outlet.amount) for outlet in created["stream"]} == {
        ("mole_flow",),
    }
    stream_flows = sorted(outlet.amount["mole_flow"]
                          for outlet in created["stream"])
    assert stream_flows == pytest.approx([2.0, 6.0])
    assert isinstance(extractor.Liquid_2, _DummyStream)
    assert isinstance(extractor.Liquid_3, _DummyStream)


def test_batch_extractor_keeps_batch_amount_semantics():
    """Batch extractor keeps inlet amount semantics as mol."""
    extractor = BatchExtractor()
    extractor.Phases = _Inlet()

    assert extractor.oper_mode == "Batch"
    assert extractor.in_flow == pytest.approx(_Inlet.moles)


@pytest.mark.parametrize("extractor_cls", [ContinuousExtractor, BatchExtractor])
def test_extractors_reject_unknown_activity_model(extractor_cls):
    """Extractor constructors reject unknown activity-model selectors."""
    with pytest.raises(ValueError, match="gamma_method must be one of"):
        extractor_cls(gamma_method="uniquac")
