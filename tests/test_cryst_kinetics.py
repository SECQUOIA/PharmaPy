import numpy as np
import pytest

from assimulo_helpers import import_module_with_assimulo_stub
from PharmaPy.Kinetics import CrystKinetics


pytestmark = pytest.mark.unit


def _primary_growth_kinetics():
    kin = CrystKinetics(
        coeff_solub=[0.1, 0.0, 0.0],
        nucl_prim=[1.0, 0.0, 1.0],
        growth=[1.0, 0.0, 1.0],
        sup_sat_type="absolute",
    )
    kin.target_idx = 0
    return kin


def _import_msmpr(monkeypatch):
    module = import_module_with_assimulo_stub(
        monkeypatch,
        "PharmaPy.Crystallizers",
        solvers={"CVode": object},
        problem={"Explicit_Problem": object},
    )
    return module.MSMPR


@pytest.mark.parametrize("conc", [0.5, np.array(0.5)])
def test_get_kinetics_accepts_scalar_target_concentration(conc):
    kin = _primary_growth_kinetics()

    nucl, growth, dissol = kin.get_kinetics(conc, 298.15, 0.5)

    assert np.ndim(nucl) == 0
    assert np.ndim(growth) == 0
    assert np.ndim(dissol) == 0
    assert nucl == pytest.approx(0.4)
    assert growth == pytest.approx(0.4)
    assert dissol == pytest.approx(0.0)


def test_get_kinetics_preserves_vector_target_selection():
    kin = _primary_growth_kinetics()
    kin.target_idx = 1

    conc = np.array([[9.0, 0.5], [8.0, 0.6]])
    temp = np.array([298.15, 298.15])
    moments = np.zeros((2, 4))

    nucl, growth, dissol = kin.get_kinetics(conc, temp, 0.5, moments)

    np.testing.assert_allclose(nucl, [0.4, 0.5])
    np.testing.assert_allclose(growth, [0.4, 0.5])
    np.testing.assert_allclose(dissol, [0.0, 0.0])


def test_get_kinetics_uses_secondary_parameters_from_vector_update():
    kin = _primary_growth_kinetics()
    kin.set_params(
        np.array([
            1.0, 0.0, 1.0,
            1.0, 0.0, 1.0, 1.0,
            1.0, 0.0, 1.0,
            0.0, 0.0, 0.0,
        ])
    )

    conc = np.array([[0.5], [0.6]])
    temp = np.array([298.15, 298.15])
    moments = np.full((2, 4), 2.0)

    prim, sec, growth, dissol = kin.get_kinetics(
        conc, temp, 0.5, moments, nucl_sec_out=True
    )

    np.testing.assert_allclose(prim, [0.4, 0.5])
    np.testing.assert_allclose(sec, [0.4, 0.5])
    np.testing.assert_allclose(growth, [0.4, 0.5])
    np.testing.assert_allclose(dissol, [0.0, 0.0])


def test_msmpr_steady_state_accepts_scalar_seed(monkeypatch):
    MSMPR = _import_msmpr(monkeypatch)

    crystallizer = MSMPR.__new__(MSMPR)
    crystallizer.vol_slurry = 1.0
    crystallizer.target_ind = 0
    crystallizer._Kinetics = _primary_growth_kinetics()

    class Solid:
        x_distrib = np.array([0.0, 1.0])
        kv = 0.0

        def getDensity(self, temp):
            return 1.0

    class Liquid:
        mass_frac = np.array([0.5])

    class Inlet:
        vol_flow = 1.0
        Liquid_1 = Liquid()

    crystallizer.Solid_1 = Solid()
    crystallizer._Inlet = Inlet()

    x_vec, f_convg, w_convg, info, final_fn = crystallizer.solve_steady_state(
        0.3, 298.15
    )

    np.testing.assert_allclose(x_vec, [0.0, 1.0])
    np.testing.assert_allclose(f_convg, [1.0, np.exp(-2.5)])
    assert w_convg == pytest.approx(0.5)
    assert info.converged
    assert final_fn == pytest.approx(0.0)
