from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("assimulo")
from PharmaPy.Drying_Model import Drying

pytestmark = pytest.mark.assimulo


def test_drying_rate_mass_basis_converts_molar_rates_with_component_mw():
    dryer = Drying(number_nodes=2, supercrit_names=["nitrogen"])
    dryer.Liquid_1 = SimpleNamespace(mw=np.array([18.0, 28.0, 46.0]))
    dry_rate = np.array([
        [2.0, 0.0, 4.0],
        [3.0, 0.0, 5.0],
    ])

    dry_rate_mass = dryer._drying_rate_mass_basis(dry_rate)

    expected = np.array([
        [0.036, 0.0, 0.184],
        [0.054, 0.0, 0.230],
    ])
    np.testing.assert_allclose(dry_rate_mass, expected)


def test_material_balance_uses_mass_drying_rate_for_saturation():
    dryer = Drying(number_nodes=2, supercrit_names=["nitrogen"])
    dryer.idx_volatiles = np.array([0, 2])
    dryer.porosity = 0.4
    dryer.rho_liq = np.array([800.0, 900.0])
    dryer.dz = np.ones(2)
    dryer.Liquid_1 = SimpleNamespace(mw=np.array([18.0, 28.0, 46.0]))

    satur = np.array([0.6, 0.8])
    temp_gas = np.array([300.0, 305.0])
    temp_sol = np.array([299.0, 304.0])
    y_gas = np.array([
        [0.02, 0.96, 0.02],
        [0.03, 0.94, 0.03],
    ])
    x_liq = np.array([
        [0.25, 0.75],
        [0.40, 0.60],
    ])
    u_gas = np.zeros(2)
    dens_gas = np.ones(2)
    dry_rate = np.array([
        [0.036, 0.0, 0.184],
        [0.054, 0.0, 0.230],
    ])
    inputs = {"mass_frac": np.array([0.01, 0.98, 0.01])}

    dsat_dt, _, _ = dryer.material_balance(
        time=0.0,
        satur=satur,
        temp_gas=temp_gas,
        temp_sol=temp_sol,
        y_gas=y_gas,
        x_liq=x_liq,
        u_gas=u_gas,
        dens_gas=dens_gas,
        dry_rate=dry_rate,
        inputs=inputs,
    )

    mass_rate = np.array([0.22, 0.284])
    expected_dsat_dt = -mass_rate / dryer.rho_liq / dryer.porosity

    np.testing.assert_allclose(dsat_dt, expected_dsat_dt)


def test_material_balance_uses_mass_drying_rate_for_gas_species():
    dryer = Drying(number_nodes=2, supercrit_names=["nitrogen"])
    dryer.idx_volatiles = np.array([0, 2])
    dryer.porosity = 0.5
    dryer.rho_liq = np.array([1000.0, 1000.0])
    dryer.dz = np.ones(2)

    satur = np.array([0.5, 0.5])
    temp_gas = np.array([300.0, 300.0])
    temp_sol = np.array([299.0, 299.0])
    y_gas = np.array([
        [0.10, 0.80, 0.10],
        [0.20, 0.70, 0.10],
    ])
    x_liq = np.array([
        [0.50, 0.50],
        [0.50, 0.50],
    ])
    u_gas = np.zeros(2)
    dens_gas = np.array([1.2, 1.5])
    dry_rate = np.array([
        [0.036, 0.0, 0.184],
        [0.054, 0.0, 0.230],
    ])
    inputs = {"mass_frac": np.array([0.05, 0.90, 0.05])}

    dsat_dt, dygas_dt, _ = dryer.material_balance(
        time=0.0,
        satur=satur,
        temp_gas=temp_gas,
        temp_sol=temp_sol,
        y_gas=y_gas,
        x_liq=x_liq,
        u_gas=u_gas,
        dens_gas=dens_gas,
        dry_rate=dry_rate,
        inputs=inputs,
    )

    epsilon_gas = dryer.porosity * (1 - satur)
    transfer_gas = (
        dry_rate / epsilon_gas[:, np.newaxis] / dens_gas[:, np.newaxis] /
        (1 - satur)[:, np.newaxis]
    )
    correction = y_gas / (1 - satur)[:, np.newaxis] * dsat_dt[:, np.newaxis]
    expected_dygas_dt = transfer_gas + correction

    np.testing.assert_allclose(dygas_dt, expected_dygas_dt)
