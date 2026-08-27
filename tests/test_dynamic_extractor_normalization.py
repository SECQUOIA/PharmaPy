"""Regression tests for DynamicExtractor mole-fraction closure."""

import numpy as np
import pytest

import PharmaPy.DynamicExtraction as dynamic_extraction


pytestmark = pytest.mark.unit


def _constant_distribution(x_light, x_heavy, temp):
    """Return asymmetric distribution coefficients for the test fixture.

    Parameters
    ----------
    x_light, x_heavy : ndarray
        Light- and heavy-phase mole fractions [-].
    temp : ndarray
        Stage temperatures [K].

    Returns
    -------
    ndarray
        Component distribution coefficients [-].
    """
    coeffs = np.array([0.4, 1.5, 1.1])  # [-]

    return np.broadcast_to(coeffs, np.shape(x_light))  # [-]


class _InletStream:
    """Minimal inlet stream exposing a molar density."""

    def __init__(self, density):
        """Store the stream density.

        Parameters
        ----------
        density : float
            Molar density [mol/L].
        """
        self.density = density  # [mol/L]

    def getDensity(self, basis="mole"):
        """Return the stream density on the requested basis.

        Parameters
        ----------
        basis : str, optional
            Density basis selector. Only ``"mole"`` is used here.

        Returns
        -------
        ndarray
            Molar density [mol/L].
        """
        return np.array(self.density)  # [mol/L]


class _Phase:
    """Minimal liquid phase exposing thermodynamic callbacks."""

    temp = 298.15  # [K]

    def getDensity(self, mole_frac=None, basis="mole", temp=None):
        """Return a constant molar density.

        Parameters
        ----------
        mole_frac : ndarray, optional
            Liquid mole fractions [-].
        basis : str, optional
            Density basis selector. Only ``"mole"`` is used here.
        temp : ndarray, optional
            Liquid temperatures [K].

        Returns
        -------
        float
            Molar density [mol/L].
        """
        return 1.0  # [mol/L]

    def getEnthalpy(self, mole_frac, temp, basis="mole"):
        """Return zero molar enthalpy for each requested stage.

        Parameters
        ----------
        mole_frac : ndarray
            Liquid mole fractions [-].
        temp : ndarray
            Liquid temperatures [K].
        basis : str, optional
            Enthalpy basis selector. Only ``"mole"`` is used here.

        Returns
        -------
        ndarray
            Molar enthalpy [J/mol].
        """
        mole_frac = np.asarray(mole_frac)  # [-]

        return np.zeros(mole_frac.shape[0])  # [J/mol]


class _BatchResult:
    """Equilibrium result fixture with normalized inlet phase fractions."""

    x_light = np.array([0.3, 0.3, 0.4])  # [-]
    x_heavy = np.array([0.2, 0.2, 0.6])  # [-]
    mol_light = 10.0  # [mol]
    mol_heavy = 8.0  # [mol]
    rho_light = 0.8  # [mol/L]
    rho_heavy = 1.2  # [mol/L]


class _BatchExtractor:
    """BatchExtractor stub returning the fixed equilibrium fixture."""

    def __init__(self, k_fun=None, gamma_method="UNIQUAC"):
        """Store constructor arguments for interface compatibility.

        Parameters
        ----------
        k_fun : callable, optional
            Equilibrium distribution callback [-].
        gamma_method : str, optional
            Activity-coefficient model selector [-].
        """
        self.k_fun = k_fun
        self.gamma_method = gamma_method
        self.result = _BatchResult()

    def solve_unit(self):
        """Skip the real batch solve.

        Returns
        -------
        None
            The prebuilt ``result`` fixture is already attached.
        """
        return None


def test_material_balances_replace_dependent_components_with_closure():
    """The dependent x and y residuals enforce mole-fraction sums."""
    extractor = dynamic_extraction.DynamicExtractor(
        num_stages=1,
        k_fun=_constant_distribution,
        eff=1.0,
    )

    x_i = np.array([[0.2, 0.3, 0.5]])  # [-]
    m_i = _constant_distribution(x_i, x_i, np.array([298.15]))  # [-]
    y_i = m_i * x_i  # [-]

    x_in = np.array([0.30, 0.30, 0.40])  # [-]
    y_in = np.array([0.20, 0.20, 0.60])  # [-]
    light_flows = np.array([10.0, 10.0])  # [mol/s]
    heavy_flows = np.array([8.0, 8.0])  # [mol/s]
    holdup_light = np.array([12.0])  # [mol]
    holdup_heavy = np.array([6.0])  # [mol]

    x_augm = np.vstack((x_in, x_i))  # [-]
    y_augm = np.vstack((y_i, y_in))  # [-]
    temp_augm = np.array([295.0, 298.15, 301.0])  # [K]
    augm_arrays = (x_augm, y_augm, temp_augm, light_flows, heavy_flows)

    residuals = extractor.material_balances(
        0.0,
        x_i=x_i,
        y_i=y_i,
        holdup_light=holdup_light,
        holdup_heavy=holdup_heavy,
        u_int=np.array([0.0]),
        temp=np.array([298.15]),
        di_sdot=None,
        augm_arrays=augm_arrays,
    )

    x_residuals, y_residuals = residuals
    component_imbalance = y_augm[1:] * heavy_flows[1:, np.newaxis] \
        + x_augm[:-1] * light_flows[:-1, np.newaxis] \
        - y_augm[:-1] * heavy_flows[:-1, np.newaxis] \
        - x_augm[1:] * light_flows[1:, np.newaxis]  # [mol/s]
    effective_holdup = holdup_light[:, np.newaxis] \
        + holdup_heavy[:, np.newaxis] * m_i  # [mol]
    expected_dxdt = component_imbalance / effective_holdup  # [1/s]

    np.testing.assert_allclose(x_residuals[:, :-1], expected_dxdt[:, :-1])
    np.testing.assert_allclose(x_residuals[:, -1], x_i.sum(axis=1) - 1.0)
    np.testing.assert_allclose(y_residuals[:, :-1],
                               m_i[:, :-1] * x_i[:, :-1] - y_i[:, :-1])
    np.testing.assert_allclose(y_residuals[:, -1], y_i.sum(axis=1) - 1.0)


def test_initialize_model_returns_closed_phase_compositions(monkeypatch):
    """The initialization correction keeps both phase fractions normalized."""
    monkeypatch.setattr(dynamic_extraction, "BatchExtractor", _BatchExtractor)

    extractor = dynamic_extraction.DynamicExtractor(
        num_stages=1,
        k_fun=_constant_distribution,
        eff=1.0,
    )
    extractor.num_comp = 3  # [-]
    extractor.Liquid_1 = _Phase()
    extractor.Inlet = {
        "feed": _InletStream(_BatchResult.rho_light),
        "solvent": _InletStream(_BatchResult.rho_heavy),
    }

    di_init = extractor.initialize_model()

    np.testing.assert_allclose(di_init["x_i"].sum(axis=1), 1.0)
    np.testing.assert_allclose(di_init["y_i"].sum(axis=1), 1.0)
