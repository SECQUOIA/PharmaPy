"""Solvent-index handling when the solvent is the first species (index 0).

`ThermoPhysicalManager.conc_to_frac`, `ThermoPhysicalManager.mass_conc_to_frac`
and `LiquidPhase.updatePhase` select their return arity from the solvent index.
A truthiness test cannot distinguish the valid index ``0`` from an absent
solvent, so these tests pin the behavior for a solvent declared as the first
species and check it against the same mixture with the species order permuted.

The self-contained three-species fixture needs no repository data. Molecular
masses and liquid densities are deliberately distinct so that an index or axis
mistake changes the numbers rather than cancelling out.
"""

import json

import numpy as np
import pytest

from PharmaPy.Phases import LiquidPhase


pytestmark = pytest.mark.unit


# Synthetic three-species liquid. Values are chosen only to keep the molar
# volumes (mw / rho_liq) well separated; they do not model a real system.
# Antoine form: log10(P/[Pa]) = A - B/(T + C), with T [K].
_SPECIES = {
    "water": {
        "mw": 18.0,  # [g/mol]
        "rho_liq": 1000.0,  # [kg/m**3]
        "t_crit": 647.0,  # [K]
        "cp_liq": [75.0],  # [J/mol/K]
        "p_vap": [8.0, 1500.0, -40.0],  # Antoine A [-], B [K], C [K]
        "delta_hvap": 40000.0,  # [J/mol]
        "tref_hvap": 350.0,  # [K]
    },
    "api": {
        "mw": 150.0,  # [g/mol]
        "rho_liq": 1200.0,  # [kg/m**3]
        "t_crit": 800.0,  # [K]
        "cp_liq": [220.0],  # [J/mol/K]
        "p_vap": [8.0, 2200.0, -40.0],  # Antoine A [-], B [K], C [K]
        "delta_hvap": 70000.0,  # [J/mol]
        "tref_hvap": 350.0,  # [K]
    },
    "impurity": {
        "mw": 60.0,  # [g/mol]
        "rho_liq": 800.0,  # [kg/m**3]
        "t_crit": 700.0,  # [K]
        "cp_liq": [130.0],  # [J/mol/K]
        "p_vap": [8.0, 1800.0, -40.0],  # Antoine A [-], B [K], C [K]
        "delta_hvap": 55000.0,  # [J/mol]
        "tref_hvap": 350.0,  # [K]
    },
}

# Solvent first: ind_solv == 0, the index a truthiness test cannot see.
SOLVENT_FIRST_ORDER = ("water", "api", "impurity")

# Same mixture with the solvent moved to a non-zero index. Species counts and
# positions differ from the order above, so a positional mistake cannot pass
# both fixtures with the same expected values.
SOLVENT_SECOND_ORDER = ("api", "water", "impurity")


def _write_thermo(tmp_path, order, filename):
    """Write a thermophysical-property file with the given species order.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory that receives the generated file.
    order : sequence of str
        Species names, in the positional order the phase should use.
    filename : str
        Name of the generated JSON file.

    Returns
    -------
    str
        Absolute path of the written file, as `LiquidPhase` expects it.
    """
    path = tmp_path / filename
    path.write_text(json.dumps({name: _SPECIES[name] for name in order}))

    return str(path)


def _molar_volumes(order):
    """Return pure-species molar volumes for a species order.

    Parameters
    ----------
    order : sequence of str
        Species names in positional order.

    Returns
    -------
    numpy.ndarray
        Pure-species molar volumes with shape ``(num_species,)`` [L/mol].

    Notes
    -----
    ``mw / rho`` converts [g/mol] over [kg/m**3] directly into [L/mol],
    because 1 kg/m**3 equals 1 g/L. This is derived from the fixture
    constants rather than from `getDensityPure`, so the expectation does not
    reuse the production conversion under test.
    """
    return np.array(
        [_SPECIES[name]["mw"] / _SPECIES[name]["rho_liq"] for name in order]
    )  # [L/mol]


def _mass_densities(order):
    """Return pure-species liquid mass densities for a species order.

    Parameters
    ----------
    order : sequence of str
        Species names in positional order.

    Returns
    -------
    numpy.ndarray
        Pure-species mass densities with shape ``(num_species,)`` [kg/m**3].
    """
    return np.array([_SPECIES[name]["rho_liq"] for name in order])  # [kg/m**3]


def _molar_masses(order):
    """Return molecular masses for a species order.

    Parameters
    ----------
    order : sequence of str
        Species names in positional order.

    Returns
    -------
    numpy.ndarray
        Molecular masses with shape ``(num_species,)`` [g/mol].
    """
    return np.array([_SPECIES[name]["mw"] for name in order])  # [g/mol]


def _expected_mole_conc(order, solvent_name, solute_conc):
    """Close the ideal-mixture volume balance on a molar basis.

    Parameters
    ----------
    order : sequence of str
        Species names in positional order.
    solvent_name : str
        Name of the solvent species.
    solute_conc : dict
        Molar concentration of each non-solvent species, keyed by name
        [mol/L].

    Returns
    -------
    numpy.ndarray
        Molar concentrations with the solvent entry filled in, shape
        ``(num_species,)`` [mol/L].

    Notes
    -----
    The solvent closes ``sum_i c_i * v_i = 1`` with ``v_i`` the pure-species
    molar volume [L/mol], so ``c_solv = (1 - sum_solutes c_i v_i) / v_solv``.
    """
    mol_vol = _molar_volumes(order)  # [L/mol]
    conc = np.array(
        [solute_conc.get(name, 0.0) for name in order]
    )  # [mol/L]

    solvent_ind = order.index(solvent_name)
    solute_volume_frac = sum(
        solute_conc[name] * mol_vol[order.index(name)] for name in solute_conc
    )  # [-]
    conc[solvent_ind] = (1.0 - solute_volume_frac) / mol_vol[solvent_ind]

    return conc


def _expected_mass_conc(order, solvent_name, solute_conc):
    """Close the ideal-mixture volume balance on a mass basis.

    Parameters
    ----------
    order : sequence of str
        Species names in positional order.
    solvent_name : str
        Name of the solvent species.
    solute_conc : dict
        Mass concentration of each non-solvent species, keyed by name
        [kg/m**3].

    Returns
    -------
    numpy.ndarray
        Mass concentrations with the solvent entry filled in, shape
        ``(num_species,)`` [kg/m**3].

    Notes
    -----
    The solvent closes ``sum_i c_i / rho_i = 1`` with ``rho_i`` the
    pure-species mass density [kg/m**3], so
    ``c_solv = (1 - sum_solutes c_i / rho_i) * rho_solv``.
    """
    rho = _mass_densities(order)  # [kg/m**3]
    conc = np.array(
        [solute_conc.get(name, 0.0) for name in order]
    )  # [kg/m**3]

    solvent_ind = order.index(solvent_name)
    solute_volume_frac = sum(
        solute_conc[name] / rho[order.index(name)] for name in solute_conc
    )  # [-]
    conc[solvent_ind] = (1.0 - solute_volume_frac) * rho[solvent_ind]

    return conc


# Solute loadings shared by the molar-basis tests, keyed by species name.
MOLE_SOLUTES = {"api": 0.4, "impurity": 0.6}  # [mol/L]

# Solute loadings shared by the mass-basis tests, keyed by species name.
MASS_SOLUTES = {"api": 60.0, "impurity": 40.0}  # [kg/m**3]


def test_conc_to_frac_returns_filled_conc_for_first_species_solvent(tmp_path):
    path = _write_thermo(tmp_path, SOLVENT_FIRST_ORDER, "solvent_first.json")
    phase = LiquidPhase(
        path, mass=1.0, mole_frac=np.array([0.98, 0.01, 0.01]),
        name_solv="water",
    )

    assert phase.ind_solv == 0

    expected_conc = _expected_mole_conc(
        SOLVENT_FIRST_ORDER, "water", MOLE_SOLUTES
    )  # [mol/L]
    expected_mole_frac = expected_conc / expected_conc.sum()  # [-]
    mw = _molar_masses(SOLVENT_FIRST_ORDER)  # [g/mol]
    expected_mass_frac = (
        expected_mole_frac * mw / np.dot(expected_mole_frac, mw)
    )  # [-]

    # Only the solute entries are read; the solvent slot is back-calculated.
    conc_in = np.array([0.0, MOLE_SOLUTES["api"], MOLE_SOLUTES["impurity"]])

    mass_frac, mole_frac, conc_out = phase.conc_to_frac(
        conc_in.copy(), solvent_ind=0
    )
    np.testing.assert_allclose(conc_out, expected_conc)
    np.testing.assert_allclose(mole_frac, expected_mole_frac)
    np.testing.assert_allclose(mass_frac, expected_mass_frac)

    mole_frac_only, conc_mole_basis = phase.conc_to_frac(
        conc_in.copy(), solvent_ind=0, basis="mole"
    )
    np.testing.assert_allclose(mole_frac_only, expected_mole_frac)
    np.testing.assert_allclose(conc_mole_basis, expected_conc)

    mass_frac_only, conc_mass_basis = phase.conc_to_frac(
        conc_in.copy(), solvent_ind=0, basis="mass"
    )
    np.testing.assert_allclose(mass_frac_only, expected_mass_frac)
    np.testing.assert_allclose(conc_mass_basis, expected_conc)


def test_mass_conc_to_frac_returns_filled_conc_for_first_species_solvent(
    tmp_path,
):
    path = _write_thermo(tmp_path, SOLVENT_FIRST_ORDER, "solvent_first.json")
    phase = LiquidPhase(
        path, mass=1.0, mole_frac=np.array([0.98, 0.01, 0.01]),
        name_solv="water",
    )

    expected_conc = _expected_mass_conc(
        SOLVENT_FIRST_ORDER, "water", MASS_SOLUTES
    )  # [kg/m**3]
    expected_mass_frac = expected_conc / expected_conc.sum()  # [-]
    mw = _molar_masses(SOLVENT_FIRST_ORDER)  # [g/mol]
    expected_mole_frac = (
        (expected_mass_frac / mw) / np.dot(expected_mass_frac, 1 / mw)
    )  # [-]

    conc_in = np.array([0.0, MASS_SOLUTES["api"], MASS_SOLUTES["impurity"]])

    mass_frac, mole_frac, conc_out = phase.mass_conc_to_frac(
        conc_in.copy(), solvent_ind=0
    )
    np.testing.assert_allclose(conc_out, expected_conc)
    np.testing.assert_allclose(mass_frac, expected_mass_frac)
    np.testing.assert_allclose(mole_frac, expected_mole_frac)

    mass_frac_only, conc_mass_basis = phase.mass_conc_to_frac(
        conc_in.copy(), solvent_ind=0, basis="mass"
    )
    np.testing.assert_allclose(mass_frac_only, expected_mass_frac)
    np.testing.assert_allclose(conc_mass_basis, expected_conc)

    mole_frac_only, conc_mole_basis = phase.mass_conc_to_frac(
        conc_in.copy(), solvent_ind=0, basis="mole"
    )
    np.testing.assert_allclose(mole_frac_only, expected_mole_frac)
    np.testing.assert_allclose(conc_mole_basis, expected_conc)


def test_liquid_phase_mole_conc_constructor_with_first_species_solvent(
    tmp_path,
):
    path = _write_thermo(tmp_path, SOLVENT_FIRST_ORDER, "solvent_first.json")

    expected_conc = _expected_mole_conc(
        SOLVENT_FIRST_ORDER, "water", MOLE_SOLUTES
    )  # [mol/L]
    expected_mole_frac = expected_conc / expected_conc.sum()  # [-]
    mw = _molar_masses(SOLVENT_FIRST_ORDER)  # [g/mol]
    expected_mass_frac = (
        expected_mole_frac * mw / np.dot(expected_mole_frac, mw)
    )  # [-]
    expected_mw_av = np.dot(expected_mole_frac, mw)  # [g/mol]

    phase = LiquidPhase(
        path,
        mass=1.0,  # [kg]
        mole_conc=np.array(
            [0.0, MOLE_SOLUTES["api"], MOLE_SOLUTES["impurity"]]
        ),
        name_solv="water",
    )

    np.testing.assert_allclose(phase.mole_conc, expected_conc)
    np.testing.assert_allclose(phase.mole_frac, expected_mole_frac)
    np.testing.assert_allclose(phase.mass_frac, expected_mass_frac)
    np.testing.assert_allclose(phase.mass_conc, expected_conc * mw)
    assert phase.mw_av == pytest.approx(expected_mw_av)


def test_liquid_phase_update_mole_conc_with_first_species_solvent(tmp_path):
    path = _write_thermo(tmp_path, SOLVENT_FIRST_ORDER, "solvent_first.json")
    phase = LiquidPhase(
        path, mass=1.0, mole_frac=np.array([0.98, 0.01, 0.01]),
        name_solv="water",
    )

    expected_conc = _expected_mole_conc(
        SOLVENT_FIRST_ORDER, "water", MOLE_SOLUTES
    )  # [mol/L]
    expected_mole_frac = expected_conc / expected_conc.sum()  # [-]
    mw = _molar_masses(SOLVENT_FIRST_ORDER)  # [g/mol]

    phase.updatePhase(
        mole_conc=np.array(
            [0.0, MOLE_SOLUTES["api"], MOLE_SOLUTES["impurity"]]
        ),
        mass=1.0,  # [kg]
    )

    np.testing.assert_allclose(phase.mole_conc, expected_conc)
    np.testing.assert_allclose(phase.mole_frac, expected_mole_frac)
    np.testing.assert_allclose(phase.mass_conc, expected_conc * mw)


def test_liquid_phase_update_mass_conc_with_first_species_solvent(tmp_path):
    path = _write_thermo(tmp_path, SOLVENT_FIRST_ORDER, "solvent_first.json")
    phase = LiquidPhase(
        path, mass=1.0, mole_frac=np.array([0.98, 0.01, 0.01]),
        name_solv="water",
    )

    expected_conc = _expected_mass_conc(
        SOLVENT_FIRST_ORDER, "water", MASS_SOLUTES
    )  # [kg/m**3]
    expected_mass_frac = expected_conc / expected_conc.sum()  # [-]
    mw = _molar_masses(SOLVENT_FIRST_ORDER)  # [g/mol]

    phase.updatePhase(
        mass_conc=np.array(
            [0.0, MASS_SOLUTES["api"], MASS_SOLUTES["impurity"]]
        ),
        mass=1.0,  # [kg]
    )

    np.testing.assert_allclose(phase.mass_conc, expected_conc)
    np.testing.assert_allclose(phase.mass_frac, expected_mass_frac)
    np.testing.assert_allclose(phase.mole_conc, expected_conc / mw)


def test_first_species_solvent_matches_reordered_species(tmp_path):
    """The same mixture must give permuted results when species are reordered.

    Reordering moves the solvent off index 0, which is the only case the
    truthiness check handled. Comparing the two orders therefore separates a
    genuine composition from a positional coincidence.
    """
    path_first = _write_thermo(
        tmp_path, SOLVENT_FIRST_ORDER, "solvent_first.json"
    )
    path_second = _write_thermo(
        tmp_path, SOLVENT_SECOND_ORDER, "solvent_second.json"
    )

    def _mole_conc_vector(order):
        return np.array(
            [MOLE_SOLUTES.get(name, 0.0) for name in order]
        )  # [mol/L]

    phase_first = LiquidPhase(
        path_first, mass=1.0, mole_conc=_mole_conc_vector(SOLVENT_FIRST_ORDER),
        name_solv="water",
    )
    phase_second = LiquidPhase(
        path_second, mass=1.0,
        mole_conc=_mole_conc_vector(SOLVENT_SECOND_ORDER),
        name_solv="water",
    )

    assert phase_first.ind_solv == 0
    assert phase_second.ind_solv == 1

    # Map the second ordering back onto the first so the same species are
    # compared, rather than the same positions.
    permutation = [SOLVENT_SECOND_ORDER.index(name)
                   for name in SOLVENT_FIRST_ORDER]

    np.testing.assert_allclose(
        phase_first.mole_conc, phase_second.mole_conc[permutation]
    )
    np.testing.assert_allclose(
        phase_first.mole_frac, phase_second.mole_frac[permutation]
    )
    np.testing.assert_allclose(
        phase_first.mass_frac, phase_second.mass_frac[permutation]
    )
    assert phase_first.mw_av == pytest.approx(phase_second.mw_av)

    # The permutation must actually reorder; otherwise the comparison above
    # would hold for the wrong reason.
    assert permutation != [0, 1, 2]


def test_absent_solvent_still_returns_fractions_only(tmp_path):
    """A phase without a declared solvent keeps the two-value return."""
    path = _write_thermo(tmp_path, SOLVENT_FIRST_ORDER, "solvent_first.json")
    phase = LiquidPhase(
        path, mass=1.0, mole_frac=np.array([0.98, 0.01, 0.01])
    )

    assert phase.ind_solv is None

    conc = np.array([50.0, 0.4, 0.6])  # [mol/L]
    mass_frac, mole_frac = phase.conc_to_frac(conc.copy())
    np.testing.assert_allclose(mole_frac, conc / conc.sum())

    mass_conc = np.array([900.0, 60.0, 40.0])  # [kg/m**3]
    mass_frac_only, mole_frac_only = phase.mass_conc_to_frac(mass_conc.copy())
    np.testing.assert_allclose(mass_frac_only, mass_conc / mass_conc.sum())

    # Single-value returns for an explicit basis must stay unwrapped.
    assert isinstance(
        phase.conc_to_frac(conc.copy(), basis="mole"), np.ndarray
    )
    assert isinstance(
        phase.mass_conc_to_frac(mass_conc.copy(), basis="mass"), np.ndarray
    )
