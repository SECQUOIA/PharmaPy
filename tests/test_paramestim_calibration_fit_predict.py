"""Regression tests for issue #78 fit and prediction paths.

Dimensional fixture names carry units, and comments call out dimensionless
intermediates as [-] where normalization or projection removes units.
"""

import numpy as np
import pytest

from PharmaPy import ParamEstim
from PharmaPy.Calibration import PCR_calibration


pytestmark = pytest.mark.unit


def test_parameter_estimation_ipopt_result_assembly_uses_base_keyword(
        monkeypatch):
    """Exercise IPOPT post-solve assembly while stubbing only the solver.

    The optional cyipopt/IPOPT dependency is replaced with a deterministic
    boundary fake, but the ParameterEstimation objective, gradient, residual,
    y_model, and covariance assembly paths remain real. Time is [s],
    concentration is [mol/L], and the fitted rate is [mol/L/s].
    """
    time_s = np.array([0.0, 1.0, 2.0])
    rate_seed_mol_l_s = 1.0
    rate_mol_l_s = 2.0
    y_model_mol_l = rate_mol_l_s * time_s
    residual_offset_mol_l = np.array([0.10, -0.05, 0.20])
    y_obs_mol_l = y_model_mol_l + residual_offset_mol_l

    def linear_model(params, x_data_s):
        """Return concentration [mol/L] from rate [mol/L/s] and time [s]."""
        return params[0] * x_data_s

    def linear_jacobian(params, x_data_s):
        """Return d(concentration)/d(rate) sensitivities with units [s]."""
        return x_data_s[np.newaxis, :]

    estimator = ParamEstim.ParameterEstimation(
        linear_model,
        param_seed=[rate_seed_mol_l_s],
        x_data=time_s,
        y_data=y_obs_mol_l,
        name_params=["rate_mol_l_s"],
        jac_fun=linear_jacobian,
    )

    def fake_minimize_ipopt(objective, params_var, jac=None, bounds=None,
                            options=None, kwargs=None):
        """Mimic IPOPT returning the solved rate [mol/L/s]."""
        optimum_mol_l_s = np.array([rate_mol_l_s])
        # Match IPOPT's solved-state callback: residuals have units of the
        # measured response before weighting, here [mol/L].
        objective(optimum_mol_l_s, **(kwargs or {}))
        return {"x": optimum_mol_l_s}

    # cyipopt/IPOPT is an optional external solver stack absent from the core
    # test lane. Patch only that boundary; objective, gradient, and covariance
    # assembly stay on the real ParameterEstimation methods.
    monkeypatch.setattr(ParamEstim, "have_cyipopt", True)
    monkeypatch.setattr(ParamEstim, "minimize_ipopt", fake_minimize_ipopt,
                        raising=False)

    opt_par_mol_l_s, covar_rate, info = estimator.optimize_fn(
        method="IPOPT", verbose=False)

    # The default identity weight matrix leaves the [mol/L] residual and [s]
    # sensitivity values numerically unchanged after sigma_inv weighting.
    expected_weighted_residuals = -residual_offset_mol_l
    expected_weighted_jacobian_s = time_s[np.newaxis, :]

    np.testing.assert_allclose(opt_par_mol_l_s, [rate_mol_l_s])
    np.testing.assert_allclose(info["fun"], expected_weighted_residuals)
    np.testing.assert_allclose(info["jac"], expected_weighted_jacobian_s)
    np.testing.assert_allclose(estimator.y_model[0].ravel(), y_model_mol_l)
    # Covariance entries correspond to rate variance units [(mol/L/s)^2].
    assert covar_rate.shape == (1, 1)


def test_pcr_predict_uses_training_centering_for_single_new_spectrum():
    """Predict a single absorbance spectrum [AU] with training statistics."""
    # Calibration predictor rows are spectra [AU] at three wavelengths.
    spectra_au = np.array([
        [0.20, 1.10, 2.40],
        [0.45, 1.35, 2.95],
        [0.80, 1.85, 3.45],
        [1.10, 2.10, 4.05],
    ])
    # Response concentrations are [g/L].
    concentration_g_l = np.array([1.2, 1.8, 2.6, 3.1])

    num_comp = 2  # number of retained principal components [-]
    calibration = PCR_calibration(spectra_au, num_comp=num_comp,
                                  standardize=True)
    calibration.get_regression(concentration_g_l, num_comp=num_comp)

    # Single prediction spectrum is in the same absorbance units [AU].
    new_spectrum_au = np.array([[0.70, 1.70, 3.20]])
    prediction_g_l = calibration.predict(new_spectrum_au)

    # Predictor centering uses the training absorbance statistics [AU], and
    # division by the training standard deviation makes the predictors [-].
    training_mean_au = spectra_au.mean(axis=0)
    training_std_au = spectra_au.std(axis=0)
    centered_new = (
        (new_spectrum_au - training_mean_au) / training_std_au
    )
    # SVD loadings and principal-component scores are dimensionless [-].
    loadings = calibration.svd_dict["V"][:, :num_comp]
    scores = centered_new @ loadings
    # Regression coefficients convert dimensionless scores [-] to [g/L], and
    # y_means is the response offset [g/L].
    regression_coeff_g_l = calibration.regression_coeff
    response_offset_g_l = calibration.y_means
    expected_g_l = scores @ regression_coeff_g_l + response_offset_g_l

    assert np.all(np.isfinite(prediction_g_l))
    np.testing.assert_allclose(prediction_g_l, expected_g_l)
