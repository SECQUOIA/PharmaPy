"""Extent-basis consistency of the reactor heat-of-reaction source term.

Issue #22: ``RxnKinetics`` normalizes the stoichiometric matrix so that the
first reactant of each reaction has coefficient 1, and the per-reaction rates
returned by ``get_rxn_rates(overall_rates=False)`` are expressed on that
normalized extent basis. The reactor energy balances instead obtain
``deltah_rxn`` from ``getHeatOfRxn`` fed the *raw* stoichiometric matrix, i.e.
per mole of reaction as written. Multiplying the two mixes bases whenever the
first reactant's coefficient differs from 1.

The regression contract exercised here is scale invariance: ``2A -> B`` with an
enthalpy of ``dh`` per mole of reaction as written describes exactly the same
chemistry as ``A -> 0.5B`` with ``dh / 2``, so both must produce the same molar
production rates *and* the same heat-release rate.

The fixture uses the pure-component database already committed for the
plug-flow reactor integration tests; species A and B have different liquid-Cp
polynomials, so the sensible-heat correction ``delta_cp`` is non-zero at the
evaluation temperature and is covered by the same invariance check.
"""

import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest


pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_INIT = REPO_ROOT / "PharmaPy" / "__init__.py"
THERMO_PATH = str(REPO_ROOT / "tests" / "integration" / "data" /
                  "pfr_test_pure_comp.json")

# Reaction written as '2A --> B' and its half-scaled equivalent 'A --> 0.5B',
# over the database species order ['A', 'B', 'C', 'solv'].
STOICH_AS_WRITTEN = np.array([[-2.0, 1.0, 0.0, 0.0]])  # [-]
STOICH_HALF_SCALED = np.array([[-1.0, 0.5, 0.0, 0.0]])  # [-]

DELTA_HRXN_AS_WRITTEN = -1.0e5  # [J/mol of reaction as written]
DELTA_HRXN_HALF_SCALED = DELTA_HRXN_AS_WRITTEN / 2  # [J/mol of reaction]

# First-order kinetics in A for both formulations, so the two objects differ
# only in how the stoichiometry (and matching enthalpy) is written.
REACTION_ORDERS = [[1.0]]  # [-]
K_PARAMS = np.array([1.0e-3])  # [1/s], pre-exponential factor
EA_PARAMS = np.array([0.0])  # [J/mol], flat temperature dependence

TREF_HRXN = 298.15  # [K]
TEMP_EVAL = 320.0  # [K], deliberately off TREF_HRXN so delta_cp is non-zero
VOL_EVAL = 2.0  # [m**3]
MOLE_CONC = np.array([1.5, 0.4, 0.0, 5.0])  # [mol/L]

STUB_HEAT_AS_WRITTEN = -100.0  # [J/mol of reaction as written]
STUB_EXTENT_RATE = 3.0  # [mol/L/s] on the normalized reaction basis
STUB_STOICH_NORMALIZATION = 2.0  # [-]
STUB_CP = np.array([100.0, 200.0])  # [J/mol/K]
STUB_MOLE_CONC = np.array([1.0, 2.0])  # [mol/L]
STUB_VOL_FLOW = 1.0  # [m**3/s]


def _stub_assimulo_modules(monkeypatch):
    """Let algebra-only reactor tests import without the optional solver."""
    assimulo = ModuleType("assimulo")

    solvers = ModuleType("assimulo.solvers")
    solvers.CVode = object
    solvers.LSODAR = object

    problem = ModuleType("assimulo.problem")
    problem.Explicit_Problem = object

    monkeypatch.setitem(sys.modules, "assimulo", assimulo)
    monkeypatch.setitem(sys.modules, "assimulo.solvers", solvers)
    monkeypatch.setitem(sys.modules, "assimulo.problem", problem)


def _prefer_source_package():
    """Avoid importing the outer checkout package named PharmaPy."""
    loaded = sys.modules.get("PharmaPy")
    loaded_path = getattr(loaded, "__file__", None)
    if loaded is not None and (
            loaded_path is None or Path(loaded_path).resolve() != PACKAGE_INIT):
        del sys.modules["PharmaPy"]

    try:
        sys.path.remove(str(REPO_ROOT))
    except ValueError:
        pass
    sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture
def pharmapy(monkeypatch):
    """Import the reactor, phase, and kinetics modules from this checkout."""
    _prefer_source_package()

    def _load():
        import PharmaPy.Kinetics as kinetics_module
        import PharmaPy.Phases as phases_module
        import PharmaPy.Reactors as reactors_module

        return kinetics_module, phases_module, reactors_module

    try:
        return _load()
    except ModuleNotFoundError as exc:
        if exc.name != "assimulo":
            raise

        _stub_assimulo_modules(monkeypatch)
        return _load()


def _build_batch_reactor(pharmapy, stoich_matrix, delta_hrxn):
    """Assemble an isothermal batch reactor for one stoichiometric writing.

    Parameters
    ----------
    pharmapy : tuple of module
        The ``(Kinetics, Phases, Reactors)`` modules from the ``pharmapy``
        fixture.
    stoich_matrix : numpy.ndarray, shape (1, 4)
        Stoichiometric coefficients [-] over the database species order.
    delta_hrxn : float
        Heat of reaction at ``TREF_HRXN`` [J/mol of reaction as written],
        consistent with ``stoich_matrix``.

    Returns
    -------
    BatchReactor
        Reactor with its liquid phase and kinetics attached and the species
        masks initialized, ready for an ``energy_balances`` call.
    """
    kinetics_module, phases_module, reactors_module = pharmapy

    liquid = phases_module.LiquidPhase(THERMO_PATH, temp=TEMP_EVAL,
                                       vol=VOL_EVAL, mole_conc=MOLE_CONC,
                                       verbose=False)

    kinetics = kinetics_module.RxnKinetics(
        path=THERMO_PATH, k_params=K_PARAMS, ea_params=EA_PARAMS,
        stoich_matrix=stoich_matrix,
        partic_species=liquid.name_species,
        params_f=REACTION_ORDERS,
        delta_hrxn=delta_hrxn, tref_hrxn=TREF_HRXN)

    reactor = reactors_module.BatchReactor(isothermal=True)
    reactor.Phases = liquid
    reactor.Kinetics = kinetics
    reactor.set_names()

    # Populated by solve_unit(); every species participates here, so the
    # inert-concentration vector is empty.
    reactor.conc_inert = liquid.mole_conc[~reactor.mask_species]  # [mol/L]

    return reactor


def _heat_generation_rate(reactor):
    """Return the reaction heat-generation rate [W] at the fixture state.

    ``energy_balances(heat_prof=True)`` reports the reaction column with the
    opposite sign of the ``source_term`` that drives ``dtemp_dt``; this helper
    undoes that so a positive value means an exothermic release.
    """
    heat_profile = reactor.energy_balances(
        0.0, MOLE_CONC, VOL_EVAL, TEMP_EVAL, TEMP_EVAL, {}, heat_prof=True)

    return -float(heat_profile[0, 0])  # [W]


class _RateBasisLiquid:
    """Minimal liquid contract for algebraic reactor energy-balance tests."""

    name_species = ["A", "B"]
    mole_conc = STUB_MOLE_CONC  # [mol/L]
    vol_flow = STUB_VOL_FLOW  # [m**3/s]

    def getCpPure(self, temp):
        """Return fixture pure-component heat capacities.

        Parameters
        ----------
        temp : float or numpy.ndarray
            Evaluation temperature [K]. The constant fixture values do not
            vary with temperature.

        Returns
        -------
        tuple
            Unused mass-basis value and molar heat capacities [J/mol/K].
        """
        return None, STUB_CP  # [J/mol/K]

    def getEnthalpy(self, temp, temp_ref, total_h=False, basis="mole"):
        """Return zero sensible enthalpy to isolate reaction heat.

        Parameters
        ----------
        temp : float or numpy.ndarray
            Evaluation temperature [K].
        temp_ref : float
            Enthalpy reference temperature [K].
        total_h : bool, optional
            Whether a total mixture enthalpy is requested. The reactor calls
            this fixture with ``False``.
        basis : str, optional
            Enthalpy basis. The reactor calls this fixture with ``"mole"``.

        Returns
        -------
        numpy.ndarray
            Zero species sensible enthalpies [J/mol].
        """
        num_temperatures = len(np.atleast_1d(temp))
        return np.zeros((num_temperatures, len(STUB_CP)))  # [J/mol]

    def getHeatOfRxn(self, stoich, temp, mask, heat_rxn_ref, tref_hrxn):
        """Return raw-basis fixture reaction enthalpy.

        Parameters
        ----------
        stoich : numpy.ndarray
            Raw stoichiometric matrix [-].
        temp : float or numpy.ndarray
            Evaluation temperature [K].
        mask : numpy.ndarray
            Participating-species mask [-].
        heat_rxn_ref : numpy.ndarray
            Reference reaction enthalpy [J/mol of reaction as written].
        tref_hrxn : float
            Reaction-enthalpy reference temperature [K].

        Returns
        -------
        numpy.ndarray
            Reaction enthalpy [J/mol of reaction as written].
        """
        num_temperatures = len(np.atleast_1d(temp))
        return np.full((num_temperatures, 1),
                       STUB_HEAT_AS_WRITTEN)  # [J/mol of reaction]


class _RateBasisKinetics:
    """Minimal normalized-rate contract for reactor energy-balance tests."""

    delta_hrxn = np.array([STUB_HEAT_AS_WRITTEN])  # [J/mol of reaction]
    stoich_matrix = np.array([[-2.0, 1.0]])  # [-]
    stoich_normalization = np.array([
        STUB_STOICH_NORMALIZATION])  # [-]
    tref_hrxn = TREF_HRXN  # [K]

    def __init__(self):
        """Initialize capture storage for the raw enthalpy handoff."""
        self.received_delta_hrxn = None  # [J/mol of reaction as written]

    def get_rxn_rates(self, conc, temp, overall_rates=False, delta_hrxn=None):
        """Return a fixed normalized extent rate [mol/L/s].

        Parameters
        ----------
        conc : numpy.ndarray
            Participating-species concentrations [mol/L].
        temp : float or numpy.ndarray
            Evaluation temperature [K].
        overall_rates : bool, optional
            Whether species rates are requested. These tests request reaction
            extent rates, so the supported value is ``False``.
        delta_hrxn : numpy.ndarray, optional
            Raw-basis reaction enthalpy [J/mol of reaction as written].

        Returns
        -------
        numpy.ndarray
            Per-reaction normalized extent rates [mol/L/s].
        """
        self.received_delta_hrxn = np.asarray(delta_hrxn)  # [J/mol of rxn]
        num_states = np.atleast_2d(conc).shape[0]
        return np.full((num_states, 1), STUB_EXTENT_RATE)  # [mol/L/s]


def _configure_rate_basis_reactor(reactor):
    """Attach deterministic liquid and kinetics stubs to a reactor.

    Parameters
    ----------
    reactor : _BaseReactor
        Reactor whose energy balance is under test.

    Returns
    -------
    _RateBasisKinetics
        Attached kinetics stub used to inspect the enthalpy handoff.
    """
    liquid = _RateBasisLiquid()
    kinetics = _RateBasisKinetics()

    reactor.Liquid_1 = liquid
    reactor._Kinetics = kinetics
    reactor._Inlet = liquid
    reactor.mask_species = np.array([True, True])

    return kinetics


def _expected_stub_source_term():
    """Return independently derived reaction heat source [W/m**3]."""
    heat_per_normalized_extent = (
        STUB_HEAT_AS_WRITTEN / STUB_STOICH_NORMALIZATION
    )  # [J/mol of normalized reaction extent]
    return -heat_per_normalized_extent * STUB_EXTENT_RATE * 1000  # [W/m**3]


def test_species_rates_are_invariant_to_stoichiometric_writing(pharmapy):
    """Guard the premise: both writings give identical molar production rates.

    If this failed, a difference in the heat term below could be attributed to
    a difference in the extent rate rather than to the mixed-basis defect.
    """
    as_written = _build_batch_reactor(pharmapy, STOICH_AS_WRITTEN,
                                      DELTA_HRXN_AS_WRITTEN)
    half_scaled = _build_batch_reactor(pharmapy, STOICH_HALF_SCALED,
                                       DELTA_HRXN_HALF_SCALED)

    rates_as_written = as_written.Kinetics.get_rxn_rates(
        MOLE_CONC, TEMP_EVAL)  # [mol/L/s]
    rates_half_scaled = half_scaled.Kinetics.get_rxn_rates(
        MOLE_CONC, TEMP_EVAL)  # [mol/L/s]

    # k * C_A with a first-order rate law and unit temperature term.
    expected_extent_rate = K_PARAMS[0] * MOLE_CONC[0]  # [mol/L/s]
    expected_rates = np.array(
        [-1.0, 0.5, 0.0, 0.0]) * expected_extent_rate  # [mol/L/s]

    np.testing.assert_allclose(rates_as_written, expected_rates, rtol=1e-10)
    np.testing.assert_allclose(rates_half_scaled, expected_rates, rtol=1e-10)


def test_heat_of_reaction_matches_normalized_extent_basis(pharmapy):
    """q_rxn must not depend on how the same reaction is written.

    ``2A --> B`` at ``dh`` and ``A --> 0.5B`` at ``dh / 2`` consume A at the
    same rate and release the same power. Issue #22: the released power scales
    with |first reactant coefficient|, so the ``2A --> B`` form reports twice
    the correct heat.
    """
    as_written = _build_batch_reactor(pharmapy, STOICH_AS_WRITTEN,
                                      DELTA_HRXN_AS_WRITTEN)
    half_scaled = _build_batch_reactor(pharmapy, STOICH_HALF_SCALED,
                                       DELTA_HRXN_HALF_SCALED)

    q_as_written = _heat_generation_rate(as_written)  # [W]
    q_half_scaled = _heat_generation_rate(half_scaled)  # [W]

    # Independent expectation: on the normalized basis the extent rate is the
    # rate of A consumption, so the released power is
    #   -dh_per_mole_A * r_A * vol, with vol in L.
    liquid = half_scaled.Liquid_1
    delta_h_per_mole_a = liquid.getHeatOfRxn(
        STOICH_HALF_SCALED, TEMP_EVAL, half_scaled.mask_species,
        np.atleast_1d(DELTA_HRXN_HALF_SCALED), TREF_HRXN)  # [J/mol of A]
    extent_rate = K_PARAMS[0] * MOLE_CONC[0]  # [mol/L/s]
    expected_q = -float(np.ravel(delta_h_per_mole_a)[0]) * extent_rate \
        * VOL_EVAL * 1000  # [W]

    assert q_half_scaled == pytest.approx(expected_q, rel=1e-10)
    assert q_as_written == pytest.approx(expected_q, rel=1e-10)


def test_batch_reported_heat_matches_temperature_source(pharmapy):
    """The scaled heat profile and temperature derivative use one source."""
    reactor = _build_batch_reactor(
        pharmapy, STOICH_AS_WRITTEN, DELTA_HRXN_AS_WRITTEN)
    reactor.heat_transfer = lambda temp, temp_ht, vol: np.zeros_like(
        np.atleast_1d(temp))  # [W]

    heat_generation = _heat_generation_rate(reactor)  # [W]
    dtemp_dt = reactor.energy_balances(
        0.0, MOLE_CONC, VOL_EVAL, TEMP_EVAL, TEMP_EVAL, {},
        heat_prof=False)  # [K/s]

    _, cp_species = reactor.Liquid_1.getCpPure(TEMP_EVAL)  # [J/mol/K]
    thermal_capacitance = VOL_EVAL * np.dot(
        MOLE_CONC * 1000, cp_species)  # [J/K]

    assert float(np.ravel(dtemp_dt)[0]) == pytest.approx(
        heat_generation / thermal_capacitance, rel=1e-10)


@pytest.mark.parametrize("reactor_name", ["CSTR", "SemibatchReactor"])
def test_tank_reactors_use_normalized_heat_with_raw_equilibrium_handoff(
        pharmapy, reactor_name):
    """CSTR and inherited semibatch paths scale only the heat source."""
    reactors_module = pharmapy[2]
    if reactor_name == "CSTR":
        reactor = reactors_module.CSTR(isothermal=True)
    else:
        reactor = reactors_module.SemibatchReactor(
            vol_tank=1.0, isothermal=True)  # [m**3]
    kinetics = _configure_rate_basis_reactor(reactor)

    reactor_volume = 2.0  # [m**3]
    inputs = {
        "Inlet": {
            "vol_flow": STUB_VOL_FLOW,
            "mole_conc": STUB_MOLE_CONC,
            "temp": TEMP_EVAL,
        }
    }
    heat_profile = reactor.energy_balances(
        0.0, STUB_MOLE_CONC, reactor_volume, TEMP_EVAL, TEMP_EVAL,
        inputs, heat_prof=True)  # [W]

    expected_source = _expected_stub_source_term() * reactor_volume  # [W]
    assert float(heat_profile[0, 0]) == pytest.approx(expected_source)
    np.testing.assert_allclose(
        kinetics.received_delta_hrxn, [[STUB_HEAT_AS_WRITTEN]])


def test_pfr_steady_energy_uses_normalized_heat_basis(pharmapy):
    """The steady PFR source uses the normalized per-reaction rate basis.

    Issue #70 owns the scalar-shape behavior of this path. The two-dimensional
    stub outputs preserve that provisional shape while independently checking
    issue #22's heat basis.
    """
    reactors_module = pharmapy[2]
    reactor = reactors_module.PlugFlowReactor(
        diam_in=1.0, num_discr=3, isothermal=False,
        adiabatic=True)  # [m], [-]
    kinetics = _configure_rate_basis_reactor(reactor)
    reactor.c_inert = np.array([])  # [mol/L]

    dtemp_dvol = reactor.energy_steady(
        STUB_MOLE_CONC, TEMP_EVAL)  # [K/m**3]

    heat_capacity_flow = STUB_VOL_FLOW * np.dot(
        STUB_CP, STUB_MOLE_CONC) * 1000  # [W/K]
    expected_dtemp_dvol = (
        _expected_stub_source_term() / heat_capacity_flow
    )  # [K/m**3]
    assert float(np.ravel(dtemp_dvol)[0]) == pytest.approx(
        expected_dtemp_dvol)
    np.testing.assert_allclose(
        kinetics.received_delta_hrxn, [[STUB_HEAT_AS_WRITTEN]])


def test_pfr_dynamic_energy_uses_normalized_heat_basis(pharmapy):
    """The dynamic finite-volume PFR source uses normalized reaction heat."""
    reactors_module = pharmapy[2]
    reactor = reactors_module.PlugFlowReactor(
        diam_in=1.0, num_discr=3, isothermal=False,
        adiabatic=True)  # [m], [-]
    _configure_rate_basis_reactor(reactor)

    reactor.vol_discr = np.array([0.0, 1.0, 2.0])  # [m**3]
    vol_diff = np.diff(reactor.vol_discr)  # [m**3]
    temp = np.full(3, TEMP_EVAL)  # [K]
    mole_conc = np.tile(STUB_MOLE_CONC, (3, 1))  # [mol/L]
    rate_i = np.full((3, 1), STUB_EXTENT_RATE)  # [mol/L/s]

    dtemp_dt = reactor.energy_balances(
        0.0, mole_conc, vol_diff, temp, STUB_VOL_FLOW, rate_i,
        heat_profile=False)  # [K/s]

    volumetric_heat_capacity = np.dot(
        STUB_CP, STUB_MOLE_CONC) * 1000  # [J/m**3/K]
    expected_dtemp_dt = (
        _expected_stub_source_term() / volumetric_heat_capacity
    )  # [K/s]
    np.testing.assert_allclose(dtemp_dt, [expected_dtemp_dt] * 2)
