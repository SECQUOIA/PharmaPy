"""Crystallization kinetics regressions for scalar target states."""

import numpy as np
import pytest

from assimulo_helpers import import_module_with_assimulo_stub
from PharmaPy.Kinetics import CrystKinetics


pytestmark = pytest.mark.unit


def _primary_growth_kinetics():
    """Create kinetics with deterministic primary nucleation and growth.

    Returns
    -------
    CrystKinetics
        Kinetics object using concentration-like inputs [kg/m**3],
        temperature [K], normalized fixture coefficients [-], primary
        nucleation [#/m**3/s], and growth [um/s].
    """
    kin = CrystKinetics(
        coeff_solub=[0.1, 0.0, 0.0],  # [kg/m**3]
        nucl_prim=[1.0, 0.0, 1.0],  # [-], normalized fixture parameters
        growth=[1.0, 0.0, 1.0],  # [-], normalized fixture parameters
        sup_sat_type="absolute",
    )
    kin.target_idx = 0  # [-]
    return kin


def _import_msmpr(monkeypatch):
    """Import MSMPR while stubbing only optional Assimulo symbols.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest cleanup fixture used by the shared optional-import helper.

    Returns
    -------
    type
        Imported ``MSMPR`` class.
    """
    module = import_module_with_assimulo_stub(
        monkeypatch,
        "PharmaPy.Crystallizers",
        solvers={"CVode": object},
        problem={"Explicit_Problem": object},
    )
    return module.MSMPR


@pytest.mark.parametrize("conc", [0.5, np.array(0.5)])
def test_get_kinetics_accepts_scalar_target_concentration(conc):
    """Scalar target concentrations return scalar kinetic rates.

    Parameters
    ----------
    conc : float or ndarray
        Target concentration fixture [kg/m**3].
    """
    kin = _primary_growth_kinetics()

    nucl, growth, dissol = kin.get_kinetics(
        conc, 298.15, 0.5)  # [#/m**3/s], [um/s], [um/s]

    assert np.ndim(nucl) == 0
    assert np.ndim(growth) == 0
    assert np.ndim(dissol) == 0
    assert nucl == pytest.approx(0.4)
    assert growth == pytest.approx(0.4)
    assert dissol == pytest.approx(0.0)


def test_get_kinetics_preserves_vector_target_selection():
    """Vector concentration inputs still select the configured target index."""
    kin = _primary_growth_kinetics()
    kin.target_idx = 1  # [-]

    conc = np.array([[9.0, 0.5], [8.0, 0.6]])  # [kg/m**3]
    temp = np.array([298.15, 298.15])  # [K]
    moments = np.zeros((2, 4))  # [-], unused normalized moment fixture

    nucl, growth, dissol = kin.get_kinetics(
        conc, temp, 0.5, moments)  # [#/m**3/s], [um/s], [um/s]

    np.testing.assert_allclose(nucl, [0.4, 0.5])
    np.testing.assert_allclose(growth, [0.4, 0.5])
    np.testing.assert_allclose(dissol, [0.0, 0.0])


def test_get_kinetics_uses_secondary_parameters_from_vector_update():
    """Vector parameter updates drive secondary nucleation rates."""
    kin = _primary_growth_kinetics()
    kin.set_params(
        np.array([
            1.0, 0.0, 1.0,
            1.0, 0.0, 1.0, 1.0,
            1.0, 0.0, 1.0,
            0.0, 0.0, 0.0,
        ])  # [-], normalized fixture parameters
    )

    conc = np.array([[0.5], [0.6]])  # [kg/m**3]
    temp = np.array([298.15, 298.15])  # [K]
    moments = np.full((2, 4), 2.0)  # [-], normalized moment fixture

    prim, sec, growth, dissol = kin.get_kinetics(
        conc, temp, 0.5, moments, nucl_sec_out=True
    )  # [#/m**3/s], [#/m**3/s], [um/s], [um/s]

    np.testing.assert_allclose(prim, [0.4, 0.5])
    np.testing.assert_allclose(sec, [0.4, 0.5])
    np.testing.assert_allclose(growth, [0.4, 0.5])
    np.testing.assert_allclose(dissol, [0.0, 0.0])


def test_msmpr_steady_state_accepts_scalar_seed(monkeypatch):
    """MSMPR steady-state solve accepts a scalar seed fraction.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest cleanup fixture used by the optional Assimulo import helper.
    """
    MSMPR = _import_msmpr(monkeypatch)

    crystallizer = MSMPR.__new__(MSMPR)
    crystallizer.vol_slurry = 1.0  # [m**3]
    crystallizer.target_ind = 0  # [-]
    crystallizer._Kinetics = _primary_growth_kinetics()

    class Solid:
        """Solid-phase fixture for MSMPR steady-state solve."""

        x_distrib = np.array([0.0, 1.0])  # [um]
        kv = 0.0  # [-]

        def getDensity(self, temp):
            """Return a constant solid density.

            Parameters
            ----------
            temp : float
                Temperature [K].

            Returns
            -------
            float
                Solid density [kg/m**3].
            """
            return 1.0  # [kg/m**3]

    class Liquid:
        """Liquid-phase fixture with one mass fraction."""

        mass_frac = np.array([0.5])  # [-]

    class Inlet:
        """Inlet fixture for the MSMPR steady-state solve."""

        vol_flow = 1.0  # [m**3/s]
        Liquid_1 = Liquid()

    crystallizer.Solid_1 = Solid()
    crystallizer._Inlet = Inlet()

    x_vec, f_convg, w_convg, info, final_fn = crystallizer.solve_steady_state(
        0.3, 298.15)  # [um], [#/m**3/um], [-], [-], [-]

    np.testing.assert_allclose(x_vec, [0.0, 1.0])
    np.testing.assert_allclose(f_convg, [1.0, np.exp(-2.5)])
    assert w_convg == pytest.approx(0.5)
    assert info.converged
    assert final_fn == pytest.approx(0.0)
