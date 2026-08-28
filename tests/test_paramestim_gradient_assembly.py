"""Regression tests for fixed-parameter and IPOPT gradient assembly.

The fixtures use concentration responses [mol/L] over time [s]. Parameter
units are stated beside each setup value because the gradients combine
state sensitivities, residual weights, and residual ordering.
"""

import numpy as np
import pytest

from PharmaPy import ParamEstim


pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    "optimize_flags, variable_params, expected_jacobian",
    [
        (
            [True, False],
            np.array([2.0]),  # [mol/L/s]
            np.array([[0.0, 1.0, 2.0]]),  # [s]
        ),
        (
            [False, True],
            np.array([0.5]),  # [mol/L]
            np.array([[1.0, 1.0, 1.0]]),  # [-]
        ),
    ],
)
def test_finite_difference_gradient_reconstructs_fixed_parameters(
        optimize_flags, variable_params, expected_jacobian):
    """Different fixed-parameter positions preserve full model parameters."""
    time_s = np.array([0.0, 1.0, 2.0])  # [s]
    param_seed = np.array([2.0, 0.5])  # [mol/L/s], [mol/L]
    observed_conc_mol_l = 2.0 * time_s + 0.5  # [mol/L]

    def affine_model(params, x_data_s):
        """Return concentration [mol/L] from rate [mol/L/s] and offset [mol/L]."""
        rate_mol_l_s, offset_mol_l = params  # [mol/L/s], [mol/L]
        return rate_mol_l_s * x_data_s + offset_mol_l

    estimator = ParamEstim.ParameterEstimation(
        affine_model,
        param_seed=param_seed,
        x_data=time_s,
        y_data=observed_conc_mol_l,
        optimize_flags=optimize_flags,
        name_params=["rate_mol_l_s", "offset_mol_l"],
    )
    estimator.optimize_flag = True
    estimator.get_objective(variable_params)

    jacobian = estimator.get_gradient(variable_params, out_array=True)

    np.testing.assert_allclose(jacobian, expected_jacobian)


def test_ipopt_gradient_uses_weighted_state_major_residuals(monkeypatch):
    """The scalar IPOPT gradient matches weighted residual objective units."""
    time_s = np.array([1.0, 2.0, 4.0])  # [s]
    rate_mol_l_s = 2.0  # [mol/L/s]
    acceleration_mol_l_s2 = 3.0  # [mol/L/s**2]
    params = np.array([rate_mol_l_s, acceleration_mol_l_s2])  # [mol/L/s], [mol/L/s**2]

    model_conc_mol_l = np.column_stack(  # [mol/L]
        (
            rate_mol_l_s * time_s,
            acceleration_mol_l_s2 * time_s**2,
        )
    )
    residuals_mol_l = np.array(  # [mol/L]
        [
            [6.0, -3.0],
            [8.0, -6.0],
            [10.0, -9.0],
        ]
    )
    observed_conc_mol_l = model_conc_mol_l - residuals_mol_l  # [mol/L]
    weight_matrix = np.diag([4.0, 9.0])  # [(mol/L)**2]

    def two_state_model(params_eval, x_data_s):
        """Return two concentration states [mol/L] at each time [s]."""
        rate_eval_mol_l_s, accel_eval_mol_l_s2 = params_eval  # [mol/L/s], [mol/L/s**2]
        return np.column_stack(
            (
                rate_eval_mol_l_s * x_data_s,
                accel_eval_mol_l_s2 * x_data_s**2,
            )
        )

    def two_state_jacobian(params_eval, x_data_s):
        """Return state-major sensitivities [s] and [s**2]."""
        sens_state_0 = np.column_stack(
            (
                x_data_s,
                np.zeros_like(x_data_s),
            )
        )  # [s], [s**2]
        sens_state_1 = np.column_stack(
            (
                np.zeros_like(x_data_s),
                x_data_s**2,
            )
        )  # [s], [s**2]
        return np.vstack((sens_state_0, sens_state_1))

    captured = {}

    def fake_minimize_ipopt(objective, params_var, jac=None, bounds=None,
                            options=None, kwargs=None):
        """Capture the IPOPT callback gradient at params_var units."""
        objective(params_var, **(kwargs or {}))
        captured["gradient"] = jac(params_var)
        return {"x": params_var}

    monkeypatch.setattr(ParamEstim, "have_cyipopt", True)
    monkeypatch.setattr(ParamEstim, "minimize_ipopt", fake_minimize_ipopt,
                        raising=False)

    estimator = ParamEstim.ParameterEstimation(
        two_state_model,
        param_seed=params,
        x_data=time_s,
        y_data=observed_conc_mol_l,
        jac_fun=two_state_jacobian,
        weight_matrix=weight_matrix,
        name_params=["rate_mol_l_s", "acceleration_mol_l_s2"],
        name_states=["linear_conc_mol_l", "quadratic_conc_mol_l"],
    )

    estimator.optimize_fn(method="IPOPT", verbose=False)

    expected_gradient = np.array([15.5, -19.0])  # [s*L/mol], [s**2*L/mol]
    np.testing.assert_allclose(captured["gradient"], expected_gradient)
