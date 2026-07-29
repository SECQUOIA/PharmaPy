"""Regression coverage for Drying condensed-phase latent heat."""

from types import SimpleNamespace

import numpy as np
import pytest

import PharmaPy.Drying_Model as drying_module
from PharmaPy.Drying_Model import Drying

pytestmark = pytest.mark.unit


def test_condensed_energy_uses_single_latent_heat_factor(monkeypatch):
    """Mass drying rates multiplied by latent heat already give full power."""
    dryer = Drying(number_nodes=4, supercrit_names=["nitrogen"])
    dryer.idx_volatiles = np.array([0, 2])  # component indices [-]
    dryer.idx_supercrit = np.array([1])  # component indices [-]
    dryer.porosity = 0.5  # [-]
    dryer.rho_sol = 1000.0  # [kg/m**3]
    dryer.cp_sol = 1000.0  # [J/kg/K]
    dryer.rho_liq = np.full(4, 800.0)  # [kg/m**3]
    dryer.h_T_j = 0.0  # [W/m**2/K]
    dryer.a_V = 1.0  # [m**2/m**3]
    dryer.h_T_loss = 0.0  # [W/m**2/K]
    dryer.cake_height = 1.0  # [m]
    dryer.T_ambient = 298.15  # [K]
    dryer.dz = np.ones(4)  # [m]

    dryer.Vapor_1 = SimpleNamespace(
        mw=np.array([28.0, 28.0, 28.0]),  # [g/mol]
        getCp=lambda temp, mass_frac, basis: np.full(4, 1000.0),  # [J/kg/K]
        getHeatVaporization=lambda temp, basis: np.array([2.0e6, 1.0e6]),  # [J/kg]
    )
    dryer.Liquid_1 = SimpleNamespace(
        getCp=lambda temp, mass_frac, basis: np.full(4, 2000.0),  # [J/kg/K]
    )

    monkeypatch.setattr(
        drying_module,
        "high_resolution_fvm",
        lambda values, boundary_cond: np.zeros(values.size + 1),  # [K]
    )

    temp_gas = np.full(4, 295.0)  # [K], matches energy_balance wet-bulb default
    temp_sol = np.array([299.0, 304.0, 309.0, 314.0])  # [K]
    satur = np.full(4, 0.5)  # [-]
    y_gas = np.array([
        [0.10, 0.80, 0.10],
        [0.20, 0.70, 0.10],
        [0.30, 0.60, 0.10],
        [0.40, 0.50, 0.10],
    ])  # gas mass fractions [-]
    x_liq = np.array([
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0, 1.0],
    ])  # liquid volatile mass fractions [-]
    dry_rate = np.array([
        [0.036, 0.0, 0.184],
        [0.054, 0.0, 0.230],
        [0.018, 0.0, 0.092],
        [0.027, 0.0, 0.046],
    ])  # [kg/m**3/s]

    _, dTcond_dt = dryer.energy_balance(
        time=0.0,
        temp_gas=temp_gas,
        temp_sol=temp_sol,
        satur=satur,
        y_gas=y_gas,
        x_liq=x_liq,
        u_gas=np.zeros(4),  # [m/s]
        rho_gas=np.ones(4),  # [kg/m**3]
        dry_rate=dry_rate,
        inputs={"temp": 295.0},  # [K]
    )

    expected_condensed_temperature_rate = np.array([
        -256000.0 / 900000.0,
        -338000.0 / 900000.0,
        -128000.0 / 900000.0,
        -100000.0 / 900000.0,
    ])  # [K/s]
    np.testing.assert_allclose(dTcond_dt, expected_condensed_temperature_rate)

    _, drying_term, _, _ = dryer.energy_balance(
        time=0.0,
        temp_gas=temp_gas,
        temp_sol=temp_sol,
        satur=satur,
        y_gas=y_gas,
        x_liq=x_liq,
        u_gas=np.zeros(4),  # [m/s]
        rho_gas=np.ones(4),  # [kg/m**3]
        dry_rate=dry_rate,
        inputs={"temp": 295.0},  # [K]
        return_terms=True,
    )
    np.testing.assert_allclose(drying_term, -expected_condensed_temperature_rate)
