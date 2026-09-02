"""Solver-free and real-backend regressions for drying-rate unit contracts."""

import numpy as np
import pytest


pytestmark = pytest.mark.unit


def test_drying_rate_mass_basis_converts_component_rates(
        drying_unit_factory):
    """Convert real-phase molar rates with each component molar mass."""
    dryer = drying_unit_factory(number_nodes=2)
    dry_rate = np.array([
        [2.0, 4.0, 0.0],
        [3.0, 5.0, 0.0],
    ])  # [mol/m**3/s]

    dry_rate_mass = dryer._drying_rate_mass_basis(dry_rate)  # [kg/m**3/s]

    expected = np.array([
        [0.036, 0.184, 0.0],
        [0.054, 0.230, 0.0],
    ])  # [kg/m**3/s], molar masses are 18, 46, and 28 g/mol
    np.testing.assert_allclose(dry_rate_mass, expected)


def test_material_balance_uses_mass_drying_rate_for_saturation(
        drying_unit_factory):
    """Use mass-basis volatile rates in the real material balance."""
    dryer = drying_unit_factory(number_nodes=2)
    dryer.idx_volatiles = np.array([0, 1])  # component indices [-]
    dryer.porosity = 0.4  # [-]
    dryer.rho_liq = np.array([800.0, 900.0])  # [kg/m**3]
    dryer.dz = np.ones(2)  # [m]

    satur = np.array([0.6, 0.8])  # [-]
    temp_gas = np.array([300.0, 305.0])  # [K]
    temp_sol = np.array([299.0, 304.0])  # [K]
    y_gas = np.array([
        [0.02, 0.02, 0.96],
        [0.03, 0.03, 0.94],
    ])  # [-]
    x_liq = np.array([
        [0.25, 0.75],
        [0.40, 0.60],
    ])  # [-]
    dry_rate = np.array([
        [0.036, 0.184, 0.0],
        [0.054, 0.230, 0.0],
    ])  # [kg/m**3/s]

    dsat_dt, _, _ = dryer.material_balance(
        time=0.0,
        satur=satur,
        temp_gas=temp_gas,
        temp_sol=temp_sol,
        y_gas=y_gas,
        x_liq=x_liq,
        u_gas=np.zeros(2),  # [m/s]
        dens_gas=np.ones(2),  # [kg/m**3]
        dry_rate=dry_rate,
        inputs={"mass_frac": np.array([0.01, 0.01, 0.98])},  # [-]
    )

    mass_rate = np.array([0.22, 0.284])  # [kg/m**3/s]
    expected_dsat_dt = -mass_rate / dryer.rho_liq / dryer.porosity  # [1/s]
    np.testing.assert_allclose(dsat_dt, expected_dsat_dt)


@pytest.mark.assimulo
def test_unit_model_hands_mass_rates_to_real_balance_path(
        drying_unit_factory):
    """The real RHS stores mass rates converted from its molar correlation."""
    pytest.importorskip("assimulo")
    dryer = drying_unit_factory(number_nodes=2)
    _, states = dryer.solve_unit(
        deltaP=5.0e4,
        runtime=1.0e-8,  # [s]
        verbose=False,
    )
    initial_state = np.asarray(states[0]).copy()  # [-] and [K]

    state_width = 3 + dryer.Liquid_1.num_species + dryer.num_volatiles  # [-]
    states_by_node = initial_state.reshape(dryer.num_nodes, state_width)
    y_gas = states_by_node[:, 1:1 + dryer.Liquid_1.num_species]  # [-]
    x_liq = states_by_node[
        :, 1 + dryer.Liquid_1.num_species:
        1 + dryer.Liquid_1.num_species + dryer.num_volatiles
    ].copy()  # [-]
    x_liq[:, -2] = 0.0  # [-], established legacy RHS state reset
    temp_cond = states_by_node[:, -1]  # [K]
    molar_rate = dryer.get_drying_rate(
        x_liq, temp_cond, y_gas, dryer.pres_gas
    )  # [mol/m**3/s]
    expected_mass_rate = dryer._drying_rate_mass_basis(
        molar_rate
    )  # [kg/m**3/s]

    model_equations = dryer.unit_model(0.0, initial_state.copy())

    assert model_equations.shape == initial_state.shape
    assert np.all(np.isfinite(model_equations))
    np.testing.assert_allclose(dryer.dry_rate, expected_mass_rate)


def test_material_balance_uses_mass_drying_rate_for_gas_species(
        drying_unit_factory):
    """Gas transfer uses one gas-holdup divisor on mass-basis rates."""
    dryer = drying_unit_factory(number_nodes=2)
    dryer.idx_volatiles = np.array([0, 1])  # component indices [-]
    dryer.porosity = 0.5  # [-]
    dryer.rho_liq = np.array([1000.0, 1000.0])  # [kg/m**3]
    dryer.dz = np.ones(2)  # [m]

    satur = np.array([0.5, 0.5])  # [-]
    y_gas = np.array([
        [0.10, 0.10, 0.80],
        [0.20, 0.10, 0.70],
    ])  # [-]
    x_liq = np.full((2, 2), 0.5)  # [-]
    dry_rate = np.array([
        [0.036, 0.184, 0.0],
        [0.054, 0.230, 0.0],
    ])  # [kg/m**3/s]
    dens_gas = np.array([1.2, 1.5])  # [kg/m**3]

    dsat_dt, dygas_dt, _ = dryer.material_balance(
        time=0.0,
        satur=satur,
        temp_gas=np.full(2, 300.0),  # [K]
        temp_sol=np.full(2, 299.0),  # [K]
        y_gas=y_gas,
        x_liq=x_liq,
        u_gas=np.zeros(2),  # [m/s]
        dens_gas=dens_gas,
        dry_rate=dry_rate,
        inputs={"mass_frac": np.array([0.01, 0.01, 0.98])},  # [-]
    )

    gas_holdup = dryer.porosity * (1.0 - satur)  # [-]
    expected_transfer = dry_rate / gas_holdup[:, None] / dens_gas[:, None]
    expected_saturation_correction = (
        y_gas / (1.0 - satur[:, None]) * dsat_dt[:, None]
    )  # [1/s]
    expected_dygas_dt = (
        expected_transfer + expected_saturation_correction
    )  # [1/s]

    np.testing.assert_allclose(dygas_dt, expected_dygas_dt)
