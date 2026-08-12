"""Drying gas mixture molar-mass regression tests.

``y_gas`` is stored on a mass basis, so the mixture molar mass is
``sum(w_i) / sum(w_i / MW_i)`` [g/mol]. For normalized mass fractions, this
reduces to ``1 / sum(w_i / MW_i)`` rather than the mass-weighted arithmetic
average ``sum(w_i * MW_i)``. The tests cover the density and energy-balance
consumers, pure-species limits, mass-/mole-basis equivalence, and invariance to
common scaling of the differential composition state.

The fixtures follow the compact synthetic convention used by the other drying
regressions. ``Vapor_1`` carries the real
``ThermoPhysicalManager.frac_to_frac`` conversion as an independent production
cross-check of the direct mass-basis calculation.
"""

import types
from types import SimpleNamespace

import numpy as np
import pytest

import PharmaPy.Drying_Model as drying_model
from PharmaPy.ThermoModule import ThermoPhysicalManager


pytestmark = pytest.mark.unit


# Universal gas constant, matching the value used by the drying model [J/mol/K].
GAS_CONSTANT = 8.314

# Species molar masses [g/mol]: water, ethanol, nitrogen. Nitrogen is the
# non-condensable drying medium; water and ethanol are the volatile pair.
MW_GAS = np.array([18.015, 46.069, 28.014])

# Node-wise gas mass fractions [-]. Rows sum to one and span a wet, water-rich
# node through a nitrogen-rich node, so the arithmetic and reciprocal-sum molar
# mass averages separate by 4-12 % across the bed.
Y_GAS = np.array([
    [0.60, 0.10, 0.30],
    [0.05, 0.15, 0.80],
    [0.25, 0.25, 0.50],
])

TEMP_GAS = np.array([313.15, 333.15, 353.15])  # [K]
TEMP_COND = np.array([310.15, 330.15, 350.15])  # [K]
PRES_GAS = np.array([101325.0, 100000.0, 98000.0])  # [Pa]
SATUR = np.array([0.45, 0.30, 0.20])  # [-]

# Liquid volatile mass fractions by node for water and ethanol [-].
X_LIQ = np.array([
    [0.30, 0.70],
    [0.40, 0.60],
    [0.50, 0.50],
])


def _mixture_molar_mass():
    """Mean molar mass of a mixture specified by mass fractions.

    Returns
    -------
    ndarray
        Node-wise mixture molar mass ``1 / sum(w_i / MW_i)`` [g/mol].
    """
    return 1 / (Y_GAS / MW_GAS).sum(axis=1)


def _expected_gas_density():
    """Ideal-gas density built from the mixture molar mass.

    Returns
    -------
    ndarray
        Node-wise gas density [kg/m**3]. The factor 1000 converts the molar
        mass from [g/mol] to [kg/mol].
    """
    return PRES_GAS * _mixture_molar_mass() / (GAS_CONSTANT * TEMP_GAS * 1000)


def _make_vapor():
    """Build the vapor collaborator shared by the regression fixtures.

    Returns
    -------
    types.SimpleNamespace
        Vapor phase exposing molar masses [g/mol], a constant mass-basis heat
        capacity [J/kg/K], a constant mass-basis heat of vaporization [J/kg],
        a constant viscosity [Pa*s], and the real ``frac_to_frac`` conversion.
    """
    vapor = SimpleNamespace(
        mw=MW_GAS,  # [g/mol]
        getViscosity=lambda temp, mass_frac: np.full(3, 1.8e-5),  # [Pa*s]
        getCp=lambda temp, mass_frac, basis: np.full(3, 1050.0),  # [J/kg/K]
        getHeatVaporization=lambda temp, basis: np.array(
            [2.40e6, 9.20e5]),  # water and ethanol [J/kg]
    )
    # The real production conversion independently cross-checks the direct
    # mass-basis calculation instead of repeating it in a test stand-in.
    vapor.frac_to_frac = types.MethodType(
        ThermoPhysicalManager.frac_to_frac, vapor)

    return vapor


def test_unit_model_gas_density_uses_mixture_molar_mass():
    """The density handed to both balances uses the mixture molar mass.

    Node 0 carries 0.60 water / 0.10 ethanol / 0.30 nitrogen by mass at
    313.15 K and 101325 Pa. Its mean molar mass is
    1 / (0.60/18.015 + 0.10/46.069 + 0.30/28.014) = 21.652 g/mol, so the density
    is 101325 * 21.652e-3 / (8.314 * 313.15) = 0.8427 kg/m**3. The mass-weighted
    arithmetic average would give 23.820 g/mol and 0.9270 kg/m**3.
    """
    dryer = drying_model.Drying(number_nodes=3, supercrit_names=["nitrogen"])

    dryer.idx_volatiles = np.array([0, 1])  # water and ethanol indices [-]
    dryer.num_volatiles = 2  # [-]
    dryer.s_inf = 0.1  # irreducible saturation [-]
    dryer.porosity = 0.5  # [-]
    dryer.rho_sol = 1500.0  # [kg/m**3]
    dryer.dPg_dz = 1.0e4  # [Pa/m]
    dryer.k_perm = 1.0e-12  # [m**2]
    dryer.pres_gas = PRES_GAS  # [Pa]
    dryer.Vapor_1 = _make_vapor()
    dryer.Liquid_1 = SimpleNamespace(
        num_species=3,
        mw=MW_GAS,  # [g/mol]
        rho_liq=np.array([997.0, 789.0, 1.14]),  # [kg/m**3]
    )

    # The drying rate is outside the density path under test; zeroing it keeps
    # the fixture free of Antoine and activity-coefficient models.
    dryer.get_drying_rate = lambda *args: np.zeros((3, 3))  # [mol/m**3/s]
    dryer.get_inputs = lambda time: {
        "Inlet": {
            "mass_frac": np.array([0.05, 0.05, 0.90]),  # [-]
            "temp": 353.15,  # [K]
        }
    }

    recorder = {}

    def material_balance(time, satur, temp_gas, temp_sol, y_gas, x_liq,
                         u_gas, dens_gas, dry_rate, inputs):
        recorder["material"] = np.array(dens_gas, copy=True)  # [kg/m**3]
        return [
            np.zeros(3),  # saturation derivative [1/s]
            np.zeros((3, 3)),  # gas mass-fraction derivative [1/s]
            np.zeros((3, 2)),  # liquid mass-fraction derivative [1/s]
        ]

    def energy_balance(time, temp_gas, temp_sol, satur, y_gas, x_liq,
                       u_gas, rho_gas, dry_rate, inputs):
        recorder["energy"] = np.array(rho_gas, copy=True)  # [kg/m**3]
        return [
            np.zeros(3),  # gas-temperature derivative [K/s]
            np.zeros(3),  # condensed-temperature derivative [K/s]
        ]

    dryer.material_balance = material_balance
    dryer.energy_balance = energy_balance

    states = np.column_stack(
        (SATUR, Y_GAS, X_LIQ, TEMP_GAS, TEMP_COND))  # [-] and [K]
    dryer.unit_model(0.0, states.ravel())

    expected = _expected_gas_density()  # [kg/m**3]
    np.testing.assert_allclose(recorder["material"], expected, rtol=1e-10)
    np.testing.assert_allclose(recorder["energy"], expected, rtol=1e-10)

    # Independent hand-computed anchor for the wet node, guarding the formula
    # above against a transcription error.
    assert recorder["material"][0] == pytest.approx(0.8427, rel=1e-3)


def test_energy_balance_specific_gas_constant_uses_mixture_molar_mass():
    """``cv = cp - R/MW`` uses the mixture molar mass, not the mass-weighted one.

    ``heat_cond`` is ``h_T_j * a_V * (Tg - Ts) / (cv * eps_gas * rho_gas)``, so
    it exposes the molar mass through ``cv``. At node 0 the mixture molar mass
    is 21.652 g/mol, giving cv = 1050 - 8.314/21.652*1000 = 666.0 J/kg/K; the
    mass-weighted average 23.820 g/mol would give 701.0 J/kg/K instead.
    """
    dryer = drying_model.Drying(number_nodes=3, supercrit_names=["nitrogen"])

    dryer.idx_volatiles = np.array([0, 1])  # water and ethanol indices [-]
    dryer.idx_supercrit = np.array([2])  # nitrogen index [-]
    dryer.porosity = 0.5  # [-]
    dryer.rho_sol = 1500.0  # [kg/m**3]
    dryer.cp_sol = 900.0  # [J/kg/K]
    dryer.rho_liq = np.array([940.0, 900.0, 880.0])  # [kg/m**3]
    dryer.h_T_j = 12.0  # [W/m**2/K]
    dryer.a_V = 2500.0  # [m**2/m**3]
    dryer.dz = np.full(3, 0.01)  # [m]
    dryer.Vapor_1 = _make_vapor()
    dryer.Liquid_1 = SimpleNamespace(
        getCp=lambda temp, mass_frac, basis: np.full(3, 2600.0),  # [J/kg/K]
    )

    rho_gas = _expected_gas_density()  # [kg/m**3], isolates the cv defect
    dry_rate = np.zeros((3, 3))  # [kg/m**3/s], removes the latent/sensible terms

    _, _, heat_cond, _ = dryer.energy_balance(
        time=0.0,
        temp_gas=TEMP_GAS,
        temp_sol=TEMP_COND,
        satur=SATUR,
        y_gas=Y_GAS,
        x_liq=X_LIQ,
        u_gas=np.zeros(3),  # [m/s]
        rho_gas=rho_gas,
        dry_rate=dry_rate,
        inputs={"temp": 353.15},  # [K]
        return_terms=True,
    )

    # Mass-basis specific gas constant R/MW; the factor 1000 converts the molar
    # mass from [g/mol] to [kg/mol].
    specific_gas_constant = GAS_CONSTANT / _mixture_molar_mass() * 1000  # [J/kg/K]
    cv_gas = 1050.0 - specific_gas_constant  # [J/kg/K]
    epsilon_gas = dryer.porosity * (1 - SATUR)  # [-]
    expected = (
        dryer.h_T_j * dryer.a_V * (TEMP_GAS - TEMP_COND)
        / (cv_gas * epsilon_gas * rho_gas)
    )  # [K/s]

    np.testing.assert_allclose(heat_cond, expected, rtol=1e-10)


@pytest.mark.parametrize("species_index", range(MW_GAS.size))
def test_gas_mixture_molar_mass_has_correct_single_species_limit(species_index):
    """A pure gas has the molar mass of its sole species [g/mol]."""
    dryer = drying_model.Drying(number_nodes=1, supercrit_names=["nitrogen"])
    dryer.Vapor_1 = _make_vapor()

    mass_fractions = np.zeros(MW_GAS.size)  # [-]
    mass_fractions[species_index] = 1.0  # [-]

    mixture_molar_mass = dryer._gas_mixture_molar_mass(
        mass_fractions
    )  # [g/mol]

    assert mixture_molar_mass == pytest.approx(MW_GAS[species_index])


def test_gas_mixture_molar_mass_matches_mole_fraction_conversion():
    """Mass- and mole-basis mixture molar masses are equivalent [g/mol]."""
    dryer = drying_model.Drying(number_nodes=3, supercrit_names=["nitrogen"])
    dryer.Vapor_1 = _make_vapor()

    mass_basis_molar_mass = dryer._gas_mixture_molar_mass(Y_GAS)  # [g/mol]
    gas_mole_fractions = dryer.Vapor_1.frac_to_frac(mass_frac=Y_GAS)  # [-]
    mole_basis_molar_mass = np.dot(gas_mole_fractions, MW_GAS)  # [g/mol]

    np.testing.assert_allclose(
        mass_basis_molar_mass,
        mole_basis_molar_mass,
        rtol=1e-12,
    )


def test_gas_mixture_molar_mass_is_invariant_to_composition_sum_drift():
    """Common scale drift does not change mixture molar mass [g/mol]."""
    dryer = drying_model.Drying(number_nodes=3, supercrit_names=["nitrogen"])
    dryer.Vapor_1 = _make_vapor()

    # Node-wise scales represent -1 %, +1 %, and +5 % drift in the freely
    # integrated gas composition sum while preserving each composition [-].
    composition_scales = np.array([0.99, 1.01, 1.05])  # [-]
    drifted_mass_fractions = (
        Y_GAS * composition_scales[:, np.newaxis]
    )  # [-]

    drifted_molar_mass = dryer._gas_mixture_molar_mass(
        drifted_mass_fractions
    )  # [g/mol]
    gas_mole_fractions = dryer.Vapor_1.frac_to_frac(
        mass_frac=drifted_mass_fractions
    )  # [-]
    mole_basis_molar_mass = np.dot(gas_mole_fractions, MW_GAS)  # [g/mol]

    np.testing.assert_allclose(
        drifted_molar_mass,
        mole_basis_molar_mass,
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        drifted_molar_mass,
        _mixture_molar_mass(),
        rtol=1e-12,
    )
