"""Drying energy-basis regression coverage.

This file uses a compact synthetic Drying fixture and calls ``unit_model`` once
instead of building a complete ``solve_unit`` drying example. That keeps the
test to a single right-hand-side evaluation while still exercising the real
``unit_model -> energy_balance`` path for issue #48. No checked-in PharmaPy
example currently constructs an end-to-end Drying run from full phase/cake
classes. External class-material examples can seed future integration coverage,
but a full Drying transient regression is deferred until the open Drying
correctness issues are resolved; this fixture documents the assumed state layout
and units for the focused energy-basis path.
"""

from types import SimpleNamespace

import numpy as np
import pytest

import PharmaPy.Drying_Model as drying_module
from PharmaPy.Drying_Model import Drying

pytestmark = pytest.mark.unit


def test_unit_model_uses_mass_drying_rate_in_energy_terms(monkeypatch):
    """Energy terms combine mass drying rates with mass-basis heat data."""
    dryer = Drying(number_nodes=3, supercrit_names=["nitrogen"])
    dryer.num_volatiles = 2  # [-]
    dryer.idx_volatiles = np.array([0, 2])  # component indices [-]
    dryer.idx_supercrit = np.array([1])  # component indices [-]
    dryer.porosity = 0.5  # [-]
    dryer.rho_sol = 1000.0  # [kg/m**3]
    dryer.cp_sol = 1000.0  # [J/kg/K]
    dryer.s_inf = 0.1  # [-]
    dryer.k_perm = 1.0  # [m**2]
    dryer.dPg_dz = 2.0  # [Pa/m]
    dryer.dz = np.ones(3)  # [m]
    dryer.pres_gas = np.array([60000.0, 122000.0, 62000.0])  # [Pa]
    dryer.CakePhase = SimpleNamespace(alpha=1.2)  # [m/kg]
    dryer.h_T_j = 0.0  # [W/m**2/K]
    dryer.a_V = 1.0  # [m**2/m**3]
    dryer.h_T_loss = 0.0  # [W/m**2/K]
    dryer.cake_height = 1.0  # [m]
    dryer.T_ambient = 298.15  # [K]

    dryer.Liquid_1 = SimpleNamespace(
        num_species=3,  # [-]
        mw=np.array([18.0, 28.0, 46.0]),  # [g/mol]
        rho_liq=np.array([800.0, 850.0, 800.0]),  # [kg/m**3]
        getCp=lambda temp, mass_frac, basis: np.array([
            2000.0,
            2000.0,
            2000.0,
        ]),  # [J/kg/K]
    )
    dryer.Vapor_1 = SimpleNamespace(
        mw=np.array([83.14, 83.14, 83.14]),  # [g/mol]
        getViscosity=lambda temp, mass_frac: np.ones(3),  # [Pa*s]
        getCp=lambda temp, mass_frac, basis: np.array([
            1000.0,
            1200.0,
            1500.0,
        ]),  # [J/kg/K]
        # One latent heat value per volatile component; the node count is
        # intentionally different so this cannot be read as per-node data.
        getHeatVaporization=lambda temp, basis: np.array([2.0e6, 1.0e6]),  # [J/kg]
    )
    dryer.get_inputs = lambda time: {
        "Inlet": {
            "mass_frac": np.array([0.05, 0.90, 0.05]),  # [-]
            "temp": 300.0,  # [K]
        }
    }

    molar_dry_rate = np.array([
        [2.0, 0.0, 4.0],
        [3.0, 0.0, 5.0],
        [1.0, 0.0, 2.0],
    ])  # [mol/m**3/s]
    monkeypatch.setattr(
        dryer,
        "get_drying_rate",
        lambda x_liq, temp_sol, y_gas, pres_gas: molar_dry_rate,
    )
    monkeypatch.setattr(
        drying_module,
        "high_resolution_fvm",
        lambda values, boundary_cond: np.zeros(values.size + 1),  # [K]
    )
    monkeypatch.setattr(
        dryer,
        "material_balance",
        lambda *args, **kwargs: [
            np.zeros(3),  # saturation derivative [1/s]
            np.zeros((3, 3)),  # gas mass-fraction derivative [1/s]
            np.zeros((3, 2)),  # liquid mass-fraction derivative [1/s]
        ],
    )

    # Columns: saturation [-], y_gas [-], x_liq [-], T_gas [K], T_cond [K].
    states = np.array([
        [0.5, 0.10, 0.80, 0.10, 0.0, 1.0, 300.0, 299.0],
        [0.5, 0.20, 0.70, 0.10, 0.0, 1.0, 305.0, 304.0],
        [0.5, 0.30, 0.60, 0.10, 0.0, 1.0, 310.0, 309.0],
    ])

    # State derivatives follow the same layout as states, with time basis [1/s].
    derivatives = dryer.unit_model(time=0.0, states=states.ravel())
    derivatives = derivatives.reshape(3, -1)

    # Sensible power [J/m**3/s] divided by gas capacity [J/m**3/K].
    # power = cp_gas [J/kg/K] * (T_gas - 295 [K]) * sum(m_dot_i) [kg/m**3/s].
    expected_gas_temperature_rate = np.array([
        1100.0 / 450.0,
        3408.0 / 1100.0,
        2475.0 / 700.0,
    ])  # [K/s]

    # The mass-rate latent powers are [256000, 338000, 128000] [J/m**3/s]. Issue #24
    # separately tracks the existing factor of two applied to those powers.
    expected_condensed_temperature_rate = np.array([
        -512000.0 / 900000.0,
        -676000.0 / 900000.0,
        -256000.0 / 900000.0,
    ])  # [K/s]

    np.testing.assert_allclose(
        derivatives[:, -2], expected_gas_temperature_rate
    )
    np.testing.assert_allclose(
        derivatives[:, -1], expected_condensed_temperature_rate
    )
