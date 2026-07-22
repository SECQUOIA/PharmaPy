"""Regression tests for raw-material accounting in ``SimulationExec``."""

import numpy as np
import pytest

from PharmaPy.SimExec import SimulationExec


pytestmark = pytest.mark.unit


class _Result:
    """Minimal solved-result test double."""

    def __init__(self, time):
        """Create a minimal result object.

        Parameters
        ----------
        time : array_like
            Simulation time grid [s].
        """
        self.time = np.asarray(time, dtype=float)  # [s]


class _DynamicInput:
    """Deterministic dynamic-input profile test double."""

    def __init__(self, **profiles):
        """Create a deterministic dynamic-input test double.

        Parameters
        ----------
        **profiles
            Dynamic profile values. Flow entries are [kg/s], [mol/s], or
            [m**3/s] according to key; temperature values are [K] and pressure
            values are [Pa].
        """
        self.profiles = profiles

    def evaluate_inputs(self, time):
        """Return dynamic raw-inlet profiles on the requested time grid.

        Parameters
        ----------
        time : ndarray
            Simulation time grid [s].

        Returns
        -------
        dict
            Flow profiles [kg/s], [mol/s], or [m**3/s] according to key;
            temperature profiles are [K] and pressure profiles are [Pa].
        """
        out = {}
        for key, value in self.profiles.items():
            if np.ndim(value) == 0:
                out[key] = np.ones_like(time, dtype=float) * value
            else:
                out[key] = np.asarray(value, dtype=float)
        return out


class _RawPhase:
    """Single-phase raw inlet test double."""

    __module__ = "PharmaPy.Phases"

    def __init__(
            self, mass_flow=0.0, mole_flow=0.0,
            mass=0.0, moles=0.0,
            mass_frac=(1.0, 0.0), mole_frac=(1.0, 0.0),
            density_mass=1000.0, density_mole=40.0):
        """Create a single-phase raw inlet.

        Parameters
        ----------
        mass_flow : float, optional
            Static mass flow [kg/s].
        mole_flow : float, optional
            Static molar flow [mol/s].
        mass : float, optional
            Batch mass or initial holdup mass [kg].
        moles : float, optional
            Batch molar amount or initial holdup amount [mol].
        mass_frac : tuple of float, optional
            Component mass fractions [-].
        mole_frac : tuple of float, optional
            Component mole fractions [-].
        density_mass : float, optional
            Mass density [kg/m**3].
        density_mole : float, optional
            Molar density [mol/L].
        """
        self.y_upstream = None
        self.DynamicInlet = None
        self.mass_flow = mass_flow  # [kg/s]
        self.mole_flow = mole_flow  # [mol/s]
        self.mass = mass  # [kg]
        self.moles = moles  # [mol]
        self.mass_frac = np.asarray(mass_frac, dtype=float)  # [-]
        self.mole_frac = np.asarray(mole_frac, dtype=float)  # [-]
        self.vol_flow = mass_flow / density_mass  # [m**3/s]
        self.vol = mass / density_mass if density_mass else 0.0  # [m**3]
        self.temp = 300.0  # [K]
        self.pres = 101325.0  # [Pa]
        self.mw_av = 50.0  # [g/mol]
        self.transferred_from_uo = False
        self._density_mass = density_mass  # [kg/m**3]
        self._density_mole = density_mole  # [mol/L]

    def getDensity(self, basis="mass", **kwargs):
        """Return the test density on the requested accounting basis.

        Parameters
        ----------
        basis : {'mass', 'mole'}, optional
            Density basis.
        **kwargs
            Ignored compatibility arguments.

        Returns
        -------
        float
            Mass density [kg/m**3] or molar density [mol/L].
        """
        if basis == "mass":
            return self._density_mass
        return self._density_mole


class _MixedRawInlet:
    """Mixed-phase raw inlet test double."""

    __module__ = "PharmaPy.MixedPhases"

    def __init__(self, phases):
        """Create a mixed inlet from phase test doubles.

        Parameters
        ----------
        phases : sequence of _RawPhase
            Component phases with flows [kg/s] and [mol/s].
        """
        self.Phases = list(phases)
        self.y_upstream = None
        self.DynamicInlet = None
        self.mass_flow = sum(phase.mass_flow for phase in phases)  # [kg/s]
        self.mole_flow = sum(phase.mole_flow for phase in phases)  # [mol/s]
        self.mass_frac = np.array([0.5, 0.5], dtype=float)  # [-]
        self.mole_frac = np.array([0.5, 0.5], dtype=float)  # [-]
        self.vol_flow = sum(phase.vol_flow for phase in phases)  # [m**3/s]
        self.temp = 300.0  # [K]
        self.pres = 101325.0  # [Pa]
        self.mw_av = 50.0  # [g/mol]
        self.transferred_from_uo = False

    def getDensity(self, basis="mass", **kwargs):
        """Return the aggregate density for compatibility.

        Parameters
        ----------
        basis : {'mass', 'mole'}, optional
            Density basis.
        **kwargs
            Ignored compatibility arguments.

        Returns
        -------
        float
            Mass density [kg/m**3] or molar density [mol/L].
        """
        return 1000.0 if basis == "mass" else 40.0


class _UnitOperation:
    """Minimal unit-operation test double."""

    __module__ = "PharmaPy.Containers"

    def __init__(self, inlet=None, time=(0.0, 10.0)):
        """Create a minimal solved unit operation.

        Parameters
        ----------
        inlet : object, optional
            Raw inlet or mixed inlet used by raw-material accounting.
        time : tuple of float, optional
            Simulated start and end times [s].
        """
        self.Inlet = inlet
        self.oper_mode = "Continuous"
        self.result = _Result(time)
        self.heat_duty = np.array([0.0, 0.0])  # [J]
        self.duty_type = np.array([0, 0], dtype=int)  # [-]
        self.outputs = {}


def _sim_with_unit(unit):
    """Build a ``SimulationExec`` instance around one unit operation.

    Parameters
    ----------
    unit : object
        Unit operation test double.

    Returns
    -------
    SimulationExec
        Simulation executor with species names [-] and one unit operation.
    """
    sim = object.__new__(SimulationExec)
    sim.NamesSpecies = ["A", "B"]
    sim.uos_instances = {"U01": unit}
    return sim


def test_get_opex_passes_raw_material_keywords_and_holdup_flag():
    """Top-level raw-material flags take precedence in OPEX accounting."""
    inlet = _RawPhase(mass_flow=2.0, mass_frac=(0.25, 0.75))
    initial_holdup = _RawPhase(mass=5.0, mass_frac=(0.4, 0.6))
    unit = _UnitOperation(inlet)
    unit.__original_phase__ = initial_holdup
    sim = _sim_with_unit(unit)

    _, raw_cost, _ = sim.GetOPEX(
        1.0,  # [USD/kg]
        include_holdups=False,
        kwargs_items={"raw_materials": {"basis": "mass"}},
    )

    expected_cost = [20.0, 5.0, 15.0]  # [USD], 2 kg/s over 10 s.

    assert "Initial_holdup" not in raw_cost.index.get_level_values(1)
    assert list(raw_cost.columns) == ["mass", "mass_A", "mass_B"]
    np.testing.assert_allclose(
        raw_cost.iloc[0][["mass", "mass_A", "mass_B"]],
        expected_cost,
    )


def test_get_opex_includes_initial_holdup_by_default():
    """OPEX includes initial holdup raw material unless disabled."""
    inlet = _RawPhase(mass_flow=2.0, mass_frac=(0.25, 0.75))
    initial_holdup = _RawPhase(mass=5.0, mass_frac=(0.4, 0.6))
    unit = _UnitOperation(inlet)
    unit.__original_phase__ = initial_holdup
    sim = _sim_with_unit(unit)

    _, raw_cost, _ = sim.GetOPEX(
        1.0,  # [USD/kg]
        kwargs_items={"raw_materials": {"basis": "mass"}},
    )

    holdup_cost = raw_cost.xs("Initial_holdup", level=1)
    expected_holdup_cost = [5.0, 2.0, 3.0]  # [USD], 5 kg at 0.4/0.6.

    assert len(holdup_cost) == 1
    np.testing.assert_allclose(
        holdup_cost.iloc[0][["mass", "mass_A", "mass_B"]],
        expected_holdup_cost,
    )


def test_mole_basis_totals_use_canonical_singular_basis():
    """Mole-basis totals use the singular ``basis='mole'`` path."""
    inlet = _RawPhase(mole_flow=4.0, mole_frac=(0.25, 0.75))
    sim = _sim_with_unit(_UnitOperation(inlet))

    raw_materials = sim.GetRawMaterials(basis="mole")
    expected_moles = [40.0, 10.0, 30.0]  # [mol], 4 mol/s over 10 s.

    assert list(raw_materials.columns) == ["moles", "moles_A", "moles_B"]
    np.testing.assert_allclose(
        raw_materials.iloc[0][["moles", "moles_A", "moles_B"]],
        expected_moles,
    )


def test_batch_raw_inlet_records_total_before_aggregation():
    """Batch raw inlets record total material before species aggregation."""
    inlet = _RawPhase(mass=6.0, mass_frac=(0.2, 0.8))
    unit = _UnitOperation(inlet)
    unit.oper_mode = "Batch"
    sim = _sim_with_unit(unit)

    raw_materials = sim.GetRawMaterials(basis="mass")
    expected_mass = [6.0, 1.2, 4.8]  # [kg], 6 kg at 0.2/0.8.

    assert list(raw_materials.columns) == ["mass", "mass_A", "mass_B"]
    np.testing.assert_allclose(
        raw_materials.iloc[0][["mass", "mass_A", "mass_B"]],
        expected_mass,
    )


def test_mixed_phase_raw_inlet_is_decomposed_by_phase():
    """Static mixed raw inlets account each phase once."""
    liquid = _RawPhase(mass_flow=1.0, mass_frac=(0.8, 0.2))
    solid = _RawPhase(mass_flow=3.0, mass_frac=(0.1, 0.9))
    mixed = _MixedRawInlet([liquid, solid])
    sim = _sim_with_unit(_UnitOperation(mixed))

    raw_materials = sim.GetRawMaterials(basis="mass", totals=False)
    expected_phase_mass = [10.0, 30.0]  # [kg], phase flows over 10 s.

    assert len(raw_materials) == 2
    np.testing.assert_allclose(
        np.sort(raw_materials["mass"].to_numpy()),
        expected_phase_mass,
    )


def test_dynamic_mixed_phase_raw_inlet_splits_total_flow_by_phase():
    """Dynamic mixed raw inlets do not integrate the total flow per phase."""
    liquid = _RawPhase(mass_flow=1.0, mass_frac=(0.8, 0.2))
    solid = _RawPhase(mass_flow=3.0, mass_frac=(0.1, 0.9))
    mixed = _MixedRawInlet([liquid, solid])
    mixed.DynamicInlet = _DynamicInput(mass_flow=4.0)  # [kg/s]
    sim = _sim_with_unit(_UnitOperation(mixed))

    raw_materials = sim.GetRawMaterials(basis="mass", totals=False)
    expected_phase_mass = [10.0, 30.0]  # [kg], 4 kg/s split 0.25/0.75.

    assert len(raw_materials) == 2
    np.testing.assert_allclose(
        np.sort(raw_materials["mass"].to_numpy()),
        expected_phase_mass,
    )


def test_get_opex_applies_nonunity_raw_cost_vector():
    """Vector raw costs price total and species columns explicitly."""
    inlet = _RawPhase(mass_flow=2.0, mass_frac=(0.25, 0.75))
    sim = _sim_with_unit(_UnitOperation(inlet))

    _, raw_cost, _ = sim.GetOPEX(
        np.array([0.0, 2.0, 3.0]),  # [USD/kg]
        include_holdups=False,
        kwargs_items={"raw_materials": {"basis": "mass"}},
    )

    expected_cost = [0.0, 10.0, 45.0]  # [USD], zero total-price column.

    np.testing.assert_allclose(
        raw_cost.iloc[0][["mass", "mass_A", "mass_B"]],
        expected_cost,
    )


def test_get_opex_rejects_raw_cost_vector_with_wrong_width():
    """Raw cost vectors must match the raw-material table width."""
    inlet = _RawPhase(mass_flow=2.0, mass_frac=(0.25, 0.75))
    sim = _sim_with_unit(_UnitOperation(inlet))

    with pytest.raises(ValueError, match="one entry per raw-material column"):
        sim.GetOPEX(
            np.array([1.0, 2.0]),  # [USD/kg]
            include_holdups=False,
            kwargs_items={"raw_materials": {"basis": "mass"}},
        )
