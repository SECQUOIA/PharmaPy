from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("assimulo")
import PharmaPy.Drying_Model as drying_module
from PharmaPy.Drying_Model import Drying

pytestmark = pytest.mark.assimulo


def test_unit_model_uses_mass_drying_rate_in_energy_terms(monkeypatch):
    """Energy terms combine mass drying rates with mass-basis heat data."""
    dryer = Drying(number_nodes=2, supercrit_names=["nitrogen"])
    dryer.num_volatiles = 2
    dryer.idx_volatiles = np.array([0, 2])
    dryer.idx_supercrit = np.array([1])
    dryer.porosity = 0.5
    dryer.rho_sol = 1000.0  # [kg/m**3]
    dryer.cp_sol = 1000.0  # [J/kg/K]
    dryer.s_inf = 0.1
    dryer.dPg_dz = 2.0
    dryer.dz = np.ones(2)
    dryer.pres_gas = np.array([60000.0, 122000.0])  # gives rho_gas [2, 4]
    dryer.CakePhase = SimpleNamespace(alpha=1.2)
    dryer.h_T_j = 0.0
    dryer.a_V = 1.0
    dryer.h_T_loss = 0.0
    dryer.cake_height = 1.0
    dryer.T_ambient = 298.15

    dryer.Liquid_1 = SimpleNamespace(
        num_species=3,
        mw=np.array([18.0, 28.0, 46.0]),  # [g/mol]
        rho_liq=np.array([800.0, 850.0, 800.0]),  # [kg/m**3]
        getCp=lambda temp, mass_frac, basis: np.array([2000.0, 2000.0]),
    )
    dryer.Vapor_1 = SimpleNamespace(
        mw=np.array([83.14, 83.14, 83.14]),  # [g/mol]
        getViscosity=lambda temp, mass_frac: np.ones(2),
        getCp=lambda temp, mass_frac, basis: np.array([1000.0, 1200.0]),
        getHeatVaporization=lambda temp, basis: np.array([2.0e6, 1.0e6]),
    )
    dryer.get_inputs = lambda time: {
        "Inlet": {
            "mass_frac": np.array([0.05, 0.90, 0.05]),
            "temp": 300.0,
        }
    }

    molar_dry_rate = np.array([
        [2.0, 0.0, 4.0],
        [3.0, 0.0, 5.0],
    ])  # [mol/m**3/s]
    monkeypatch.setattr(
        dryer,
        "get_drying_rate",
        lambda x_liq, temp_sol, y_gas, pres_gas: molar_dry_rate,
    )
    monkeypatch.setattr(
        drying_module,
        "high_resolution_fvm",
        lambda values, boundary_cond: np.zeros(values.size + 1),
    )
    monkeypatch.setattr(
        dryer,
        "material_balance",
        lambda *args, **kwargs: [
            np.zeros(2),
            np.zeros((2, 3)),
            np.zeros((2, 2)),
        ],
    )

    states = np.array([
        [0.5, 0.10, 0.80, 0.10, 0.0, 1.0, 300.0, 299.0],
        [0.5, 0.20, 0.70, 0.10, 0.0, 1.0, 305.0, 304.0],
    ])

    derivatives = dryer.unit_model(time=0.0, states=states.ravel())
    derivatives = derivatives.reshape(2, -1)

    # Sensible power is cp_gas * (T_gas - 295 K) * sum(m_dot_i).
    expected_gas_temperature_rate = np.array([1100.0 / 450.0, 3408.0 / 1100.0])

    # The mass-rate latent powers are [256000, 338000] J/m**3/s. Issue #24
    # separately tracks the existing factor of two applied to those powers.
    expected_condensed_temperature_rate = np.array([
        -512000.0 / 900000.0,
        -676000.0 / 900000.0,
    ])

    np.testing.assert_allclose(
        derivatives[:, -2], expected_gas_temperature_rate
    )
    np.testing.assert_allclose(
        derivatives[:, -1], expected_condensed_temperature_rate
    )
