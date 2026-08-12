"""Shape contract of ``VaporPhase.getHeatVaporization`` and its Drying consumer.

Regression coverage for issue #101. The mass-basis branch used to reassign
``deltahvap`` to its subcritical column subset, so a phase holding any species
that is supercritical at every requested temperature returned fewer columns on
a mass basis than on a mole basis. Mixture weighting in
``VaporPhase.getEnthalpy`` then failed, and ``Drying.energy_balance`` only
worked while that accidental narrowing happened to match its volatile count.

Expected latent heats are derived here with a scalar loop over one temperature
at a time, independent of the vectorised two-dimensional fancy indexing under
test, which is where the defect lived. The Drying fixture is synthetic and
calls ``energy_balance`` directly, following the sibling drying tests; its
carrier latent heat carries a deliberately non-zero sentinel so that any
consumption of the non-condensable column is numerically obvious.
"""

import json
from types import SimpleNamespace

import numpy as np
import pytest

import PharmaPy.Drying_Model as drying_module
from PharmaPy.Drying_Model import Drying
from PharmaPy.Phases import VaporPhase

pytestmark = pytest.mark.unit


# Two species whose critical temperatures straddle the probe temperatures, so
# one species is supercritical while the other still condenses.
THERMO_TWO_SPECIES = {
    "light": {
        "mw": 18.0,
        "t_crit": 650.0,
        "rho_liq": 1000.0,
        "cp_liq": [75.0],
        "cp_vapor": [35.0],
        "p_vap": [8.0, 1500.0, -40.0],
        "delta_hvap": 40000.0,
        "tref_hvap": 350.0,
    },
    "heavy": {
        "mw": 100.0,
        "t_crit": 700.0,
        "rho_liq": 900.0,
        "cp_liq": [150.0],
        "cp_vapor": [80.0],
        "p_vap": [8.0, 1800.0, -40.0],
        "delta_hvap": 60000.0,
        "tref_hvap": 350.0,
    },
}

MOLAR_MASS = np.array([18.0, 100.0])  # [g/mol]
T_CRIT = np.array([650.0, 700.0])  # [K]
DELTA_HVAP_REF = np.array([40000.0, 60000.0])  # [J/mol] at TREF_HVAP
TREF_HVAP = np.array([350.0, 350.0])  # [K]

# Exponent of the Watson correlation used by VaporPhase.getHeatVaporization to
# extrapolate the tabulated latent heat away from its reference temperature.
WATSON_EXPONENT = 0.38  # [-]

G_PER_KG = 1000.0  # [g/kg], molar mass is tabulated in g/mol


def _latent_heat_mole_reference(temp_scalar):
    """Latent heat of every species at one temperature, on a molar basis.

    Scalar reference implementation of the Watson correlation, written as a
    loop over species so that it shares no indexing logic with the vectorised
    implementation under test.

    Parameters
    ----------
    temp_scalar : float
        Temperature at which the latent heat is evaluated [K].

    Returns
    -------
    ndarray
        Latent heat of each species [J/mol]. Entries for species that are
        supercritical at ``temp_scalar`` are zero, because a supercritical
        species cannot condense.
    """
    latent_heat = np.zeros(len(MOLAR_MASS))  # [J/mol]
    for species in range(len(MOLAR_MASS)):
        if temp_scalar < T_CRIT[species]:
            reduced = ((T_CRIT[species] - temp_scalar)
                       / (T_CRIT[species] - TREF_HVAP[species]))  # [-]
            latent_heat[species] = (DELTA_HVAP_REF[species]
                                    * reduced**WATSON_EXPONENT)

    return latent_heat


@pytest.fixture
def vapor_phase(tmp_path):
    """Two-species ``VaporPhase`` built from the straddling-t_crit fixture."""
    path = tmp_path / "thermo_two_species.json"
    path.write_text(json.dumps(THERMO_TWO_SPECIES))

    return VaporPhase(str(path), temp=350.0, moles=1.0,
                      mole_frac=np.array([0.5, 0.5]), check_input=False)


def test_multitemperature_mass_basis_keeps_supercritical_columns(vapor_phase):
    """A supercritical species keeps a zero column instead of being dropped.

    Three temperatures against two species keep the temperature and species
    axes different lengths, so a transposed or mis-sliced result cannot match
    by coincidence.
    """
    # All three lie above t_crit of "light" and below t_crit of "heavy".
    temps = np.array([660.0, 670.0, 680.0])  # [K]

    latent_mass = vapor_phase.getHeatVaporization(temps, basis="mass")
    latent_mole = vapor_phase.getHeatVaporization(temps, basis="mole")

    assert latent_mass.shape == (len(temps), len(MOLAR_MASS))
    assert latent_mass.shape == latent_mole.shape

    expected_mole = np.vstack([_latent_heat_mole_reference(t) for t in temps])
    expected_mass = expected_mole / MOLAR_MASS * G_PER_KG  # [J/kg]

    np.testing.assert_allclose(latent_mole, expected_mole, rtol=1e-12)
    np.testing.assert_allclose(latent_mass, expected_mass, rtol=1e-12)

    # "light" is supercritical at both temperatures and cannot condense.
    np.testing.assert_array_equal(latent_mass[:, 0], np.zeros(len(temps)))
    assert np.all(latent_mass[:, 1] > 0.0)


def test_single_temperature_mass_basis_keeps_supercritical_columns(vapor_phase):
    """The scalar-temperature branch preserves the species axis as well."""
    temp = 660.0  # [K], above t_crit of "light" only

    latent_mass = vapor_phase.getHeatVaporization(temp, basis="mass")

    assert latent_mass.shape == (len(MOLAR_MASS),)

    expected_mass = (_latent_heat_mole_reference(temp)
                     / MOLAR_MASS * G_PER_KG)  # [J/kg]
    np.testing.assert_allclose(latent_mass, expected_mass, rtol=1e-12)
    assert latent_mass[0] == 0.0


def test_getenthalpy_mass_basis_weights_full_species_axis(vapor_phase):
    """Mixture weighting no longer fails when a species is supercritical.

    This is the failure reported in issue #101: the mass-basis latent heat
    reached ``np.dot`` with fewer columns than the mass-fraction vector had
    entries. The temperature count has to equal the species count here, since
    ``getEnthalpy`` separately compares ``temp`` against ``t_crit``
    elementwise and cannot yet accept a differing number of temperatures.
    """
    temps = np.array([660.0, 680.0])  # [K]

    enthalpy_mass = vapor_phase.getEnthalpy(temp=temps, basis="mass")

    assert enthalpy_mass.shape == (len(temps),)
    assert np.all(np.isfinite(enthalpy_mass))


@pytest.mark.parametrize("temps", [
    np.array([400.0, 420.0, 440.0]),  # [K], every species subcritical
    np.array([400.0]),  # [K], single-element array
])
def test_all_subcritical_mass_basis_is_unchanged(vapor_phase, temps):
    """The previously working all-subcritical path keeps its exact values.

    With no supercritical species the old expression indexed every column, so
    the in-place conversion must reproduce it exactly rather than merely
    closely.
    """
    latent_mass = vapor_phase.getHeatVaporization(temps, basis="mass")
    latent_mole = vapor_phase.getHeatVaporization(temps, basis="mole")

    # Pre-fix expression, evaluated on the same molar intermediate.
    every_species = np.arange(len(MOLAR_MASS))
    previous = (np.atleast_2d(latent_mole)[:, every_species]
                / MOLAR_MASS[every_species] * G_PER_KG)  # [J/kg]

    np.testing.assert_array_equal(np.atleast_2d(latent_mass), previous)
    assert np.all(latent_mass > 0.0)


def test_energy_balance_uses_only_volatile_latent_heat_columns(monkeypatch):
    """Drying pairs each volatile's drying rate with that species' latent heat.

    The fixture keeps the node count, species count, and volatile count all
    different (4, 3, 2) so that an axis or ordering mistake cannot pass by
    coincidence, and gives the non-condensable carrier a large sentinel latent
    heat that the real correlation would report as zero.
    """
    num_nodes = 4  # [-]
    dryer = Drying(number_nodes=num_nodes, supercrit_names=["nitrogen"])
    dryer.idx_volatiles = np.array([0, 2])  # species indices [-]
    dryer.idx_supercrit = np.array([1])  # species indices [-]
    dryer.porosity = 0.5  # [-]
    dryer.rho_sol = 1000.0  # [kg/m**3]
    dryer.cp_sol = 1000.0  # [J/kg/K]
    dryer.rho_liq = np.full(num_nodes, 800.0)  # [kg/m**3]
    dryer.h_T_j = 0.0  # [W/m**2/K], isolates the latent term
    dryer.a_V = 1.0  # [m**2/m**3]
    dryer.h_T_loss = 0.0  # [W/m**2/K]
    dryer.cake_height = 1.0  # [m]
    dryer.T_ambient = 298.15  # [K]
    dryer.dz = np.ones(num_nodes)  # [m]

    volatile_latent_heat = np.array([2.0e6, 1.0e6])  # [J/kg], species 0 and 2
    # Species 1 is the non-condensable carrier. The real correlation returns
    # zero for it; this sentinel makes any use of that column obvious.
    carrier_sentinel = 7.0e6  # [J/kg]
    observed_temps = []

    def latent_heat_stub(temp, basis):
        """Full-width latent heat, one row per node and one column per species."""
        observed_temps.append(np.atleast_1d(temp))
        row = np.array([volatile_latent_heat[0], carrier_sentinel,
                        volatile_latent_heat[1]])  # [J/kg]
        return np.tile(row, (len(np.atleast_1d(temp)), 1))

    dryer.Vapor_1 = SimpleNamespace(
        mw=np.array([28.0, 28.0, 28.0]),  # [g/mol]
        getCp=lambda temp, mass_frac, basis: np.full(num_nodes, 1000.0),
        getHeatVaporization=latent_heat_stub,
    )
    dryer.Liquid_1 = SimpleNamespace(
        getCp=lambda temp, mass_frac, basis: np.full(num_nodes, 2000.0),
    )

    monkeypatch.setattr(
        drying_module,
        "high_resolution_fvm",
        lambda values, boundary_cond: np.zeros(values.size + 1),  # [K]
    )

    temp_sol = np.array([299.0, 304.0, 309.0, 314.0])  # [K]
    dry_rate = np.array([
        [0.036, 0.5, 0.184],
        [0.054, 0.5, 0.230],
        [0.018, 0.5, 0.092],
        [0.027, 0.5, 0.046],
    ])  # [kg/m**3/s]; the carrier column is non-zero on purpose

    _, dTcond_dt = dryer.energy_balance(
        time=0.0,
        temp_gas=np.full(num_nodes, 295.0),  # [K]
        temp_sol=temp_sol,
        satur=np.full(num_nodes, 0.5),  # [-]
        y_gas=np.tile(np.array([0.10, 0.80, 0.10]), (num_nodes, 1)),  # [-]
        x_liq=np.tile(np.array([0.0, 1.0]), (num_nodes, 1)),  # [-]
        u_gas=np.zeros(num_nodes),  # [m/s]
        rho_gas=np.ones(num_nodes),  # [kg/m**3]
        dry_rate=dry_rate,
        inputs={"temp": 295.0},  # [K]
    )

    # The latent heat is requested once per node, at the condensed-phase
    # temperatures, not once per volatile species.
    assert observed_temps
    np.testing.assert_allclose(observed_temps[0], temp_sol)

    # Latent power sums only the volatile species [J/m**3/s].
    expected_latent_power = (dry_rate[:, [0, 2]]
                             * volatile_latent_heat).sum(axis=1)
    # Solid and liquid heat capacity per unit cake volume [J/m**3/K].
    heat_capacity = (1000.0 * (1 - 0.5) * 1000.0
                     + 0.5 * 0.5 * 2000.0 * 800.0)
    expected_rate = -expected_latent_power / heat_capacity  # [K/s]

    np.testing.assert_allclose(dTcond_dt, expected_rate)

    # A run that consumed the carrier column would be far larger in magnitude.
    leaked = -(dry_rate * np.array([volatile_latent_heat[0], carrier_sentinel,
                                    volatile_latent_heat[1]])
               ).sum(axis=1) / heat_capacity  # [K/s]
    assert np.all(np.abs(dTcond_dt) < np.abs(leaked) / 2)
