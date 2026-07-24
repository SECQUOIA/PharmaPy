import numpy as np
import pytest

from PharmaPy import ParamEstim
from PharmaPy.Calibration import PCR_calibration


pytestmark = pytest.mark.unit


def test_parameter_estimation_ipopt_result_assembly_uses_base_keyword(
        monkeypatch):
    time_s = np.array([0.0, 1.0, 2.0])
    rate_mol_l_s = 2.0
    y_obs_mol_l = rate_mol_l_s * time_s + np.array([0.10, -0.05, 0.20])

    def linear_model(params, x_data_s):
        return params[0] * x_data_s

    estimator = ParamEstim.ParameterEstimation(
        linear_model,
        param_seed=[1.0],
        x_data=time_s,
        y_data=y_obs_mol_l,
        name_params=["rate_mol_l_s"],
    )

    def fake_minimize_ipopt(objective, params_var, jac=None, bounds=None,
                            options=None, kwargs=None):
        optimum = np.array([rate_mol_l_s])
        # Match IPOPT's solved-state callback: residuals have units of the
        # measured response before weighting, here [mol/L].
        objective(optimum, **(kwargs or {}))
        return {"x": optimum}

    def final_gradient(params, out_array=False, set_self=True):
        # d(y_model)/d(rate) has units [s], so weighted residual derivatives
        # preserve the expected parameter-estimation covariance dimensions.
        jacobian = time_s[np.newaxis, :]
        if out_array:
            return jacobian
        return np.zeros(1)

    monkeypatch.setattr(ParamEstim, "have_cyipopt", True)
    monkeypatch.setattr(ParamEstim, "minimize_ipopt", fake_minimize_ipopt,
                        raising=False)
    monkeypatch.setattr(estimator, "get_gradient", final_gradient)

    opt_par, covar, info = estimator.optimize_fn(method="IPOPT",
                                                 verbose=False)

    np.testing.assert_allclose(opt_par, [rate_mol_l_s])
    np.testing.assert_allclose(info["fun"], [-0.10, 0.05, -0.20])
    np.testing.assert_allclose(estimator.y_model[0].ravel(),
                               rate_mol_l_s * time_s)
    assert covar.shape == (1, 1)


def test_pcr_predict_uses_training_centering_for_single_new_spectrum():
    spectra_au = np.array([
        [0.20, 1.10, 2.40],
        [0.45, 1.35, 2.95],
        [0.80, 1.85, 3.45],
        [1.10, 2.10, 4.05],
    ])
    concentration_g_l = np.array([1.2, 1.8, 2.6, 3.1])

    calibration = PCR_calibration(spectra_au, num_comp=2, standardize=True)
    calibration.get_regression(concentration_g_l, num_comp=2)

    new_spectrum_au = np.array([[0.70, 1.70, 3.20]])
    prediction_g_l = calibration.predict(new_spectrum_au)

    # Predictor centering uses the training absorbance statistics [AU].
    centered_new = (
        (new_spectrum_au - spectra_au.mean(axis=0)) / spectra_au.std(axis=0)
    )
    # Standardized predictors and principal-component scores are
    # dimensionless [-]; regression coefficients restore [g/L].
    scores = centered_new @ calibration.svd_dict["V"][:, :2]
    expected_g_l = scores @ calibration.regression_coeff + calibration.y_means

    assert np.all(np.isfinite(prediction_g_l))
    np.testing.assert_allclose(prediction_g_l, expected_g_l)
