"""Regression tests for multistage DynamicExtractor efficiency corrections."""

import numpy as np
import pytest

import PharmaPy.DynamicExtraction as dynamic_extraction


pytestmark = pytest.mark.unit


_K_STAGE = np.array([
    [0.60, 1.20, 1.50],
    [1.80, 0.75, 1.05],
])  # [-]


def _stage_distribution(x_light, x_heavy, temp):
    """Return stage-specific distribution coefficients.

    Parameters
    ----------
    x_light, x_heavy : ndarray
        Light- and heavy-phase mole fractions [-].
    temp : ndarray
        Stage temperatures [K].

    Returns
    -------
    ndarray
        Stage-specific distribution coefficients ``K_i`` [-].
    """
    return np.broadcast_to(_K_STAGE, np.shape(x_light))  # [-]


def test_material_balances_apply_stage_aligned_efficiency_correction():
    """Later-stage correction uses same-stage ``K_i`` and heavy holdup."""
    extractor = dynamic_extraction.DynamicExtractor(
        num_stages=2,
        k_fun=_stage_distribution,
        eff=0.75,
    )

    x_i = np.array([
        [0.20, 0.35, 0.45],
        [0.32, 0.28, 0.40],
    ])  # [-]
    y_i = np.array([
        [0.25, 0.30, 0.45],
        [0.18, 0.42, 0.40],
    ])  # [-]
    x_in = np.array([0.30, 0.30, 0.40])  # [-]
    y_in = np.array([0.22, 0.38, 0.40])  # [-]
    x_augm = np.vstack((x_in, x_i))  # [-]
    y_augm = np.vstack((y_i, y_in))  # [-]
    temp_augm = np.array([296.15, 298.15, 299.15, 301.15])  # [K]
    light_flows = np.array([5.0, 4.0, 3.5])  # [mol/s]
    heavy_flows = np.array([2.5, 3.0, 3.2])  # [mol/s]
    holdup_light = np.array([10.0, 14.0])  # [mol]
    holdup_heavy = np.array([4.0, 9.0])  # [mol]
    augm_arrays = (x_augm, y_augm, temp_augm, light_flows, heavy_flows)

    material_residuals, equilibrium_residuals = extractor.material_balances(
        0.0,
        x_i=x_i,
        y_i=y_i,
        holdup_light=holdup_light,
        holdup_heavy=holdup_heavy,
        u_int=np.array([0.0, 0.0]),  # [J]
        temp=np.array([298.15, 299.15]),  # [K]
        di_sdot=None,
        augm_arrays=augm_arrays,
    )

    expected_material = np.array([
        [0.0465909090909091, 0.03719512195121951, 0.0],
        [0.00268513789581205, 0.01998647932131495, 0.0],
    ])  # first two columns [1/s], final column [-]
    expected_equilibrium = np.array([
        [-0.1500, 0.1400, 0.0],
        [0.4680, -0.2275, 0.0],
    ])  # [-]
    # The second-stage rates include
    # H_E,2 * (K_i,2 / eff) * (1 - eff) / div_2 * dx_i,1/dt.
    # The unequal K_i and H_E values guard against using stage 1 quantities.

    np.testing.assert_allclose(material_residuals, expected_material)
    np.testing.assert_allclose(equilibrium_residuals, expected_equilibrium)
