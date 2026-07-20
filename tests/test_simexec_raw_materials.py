import numpy as np
import pytest

from PharmaPy.SimExec import SimulationExec


pytestmark = pytest.mark.unit


class _Result:
    def __init__(self, time):
        self.time = np.asarray(time, dtype=float)


class _RawPhase:
    __module__ = "PharmaPy.Phases"

    def __init__(
            self, mass_flow=0.0, mole_flow=0.0,
            mass=0.0, moles=0.0,
            mass_frac=(1.0, 0.0), mole_frac=(1.0, 0.0),
            density_mass=1000.0, density_mole=40.0):
        self.y_upstream = None
        self.DynamicInlet = None
        self.mass_flow = mass_flow  # kg/s
        self.mole_flow = mole_flow  # mol/s
        self.mass = mass  # kg
        self.moles = moles  # mol
        self.mass_frac = np.asarray(mass_frac, dtype=float)  # [-]
        self.mole_frac = np.asarray(mole_frac, dtype=float)  # [-]
        self.vol_flow = mass_flow / density_mass  # m**3/s
        self.vol = mass / density_mass if density_mass else 0.0  # m**3
        self.temp = 300.0  # K
        self.pres = 101325.0  # Pa
        self.mw_av = 50.0  # g/mol
        self.transferred_from_uo = False
        self._density_mass = density_mass  # kg/m**3
        self._density_mole = density_mole  # mol/L

    def getDensity(self, basis="mass", **kwargs):
        if basis == "mass":
            return self._density_mass
        return self._density_mole


class _MixedRawInlet:
    __module__ = "PharmaPy.MixedPhases"

    def __init__(self, phases):
        self.Phases = list(phases)
        self.y_upstream = None
        self.DynamicInlet = None
        self.mass_flow = sum(phase.mass_flow for phase in phases)  # kg/s
        self.mole_flow = sum(phase.mole_flow for phase in phases)  # mol/s
        self.mass_frac = np.array([0.5, 0.5], dtype=float)  # [-]
        self.mole_frac = np.array([0.5, 0.5], dtype=float)  # [-]
        self.vol_flow = sum(phase.vol_flow for phase in phases)  # m**3/s
        self.temp = 300.0  # K
        self.pres = 101325.0  # Pa
        self.mw_av = 50.0  # g/mol
        self.transferred_from_uo = False

    def getDensity(self, basis="mass", **kwargs):
        return 1000.0 if basis == "mass" else 40.0


class _UnitOperation:
    __module__ = "PharmaPy.Containers"

    def __init__(self, inlet=None, time=(0.0, 10.0)):
        self.Inlet = inlet
        self.oper_mode = "Continuous"
        self.result = _Result(time)
        self.heat_duty = np.array([0.0, 0.0])  # J
        self.duty_type = np.array([0, 0], dtype=int)
        self.outputs = {}


def _sim_with_unit(unit):
    sim = object.__new__(SimulationExec)
    sim.NamesSpecies = ["A", "B"]
    sim.uos_instances = {"U01": unit}
    return sim


def test_get_opex_passes_raw_material_keywords_and_holdup_flag():
    inlet = _RawPhase(mass_flow=2.0, mass_frac=(0.25, 0.75))
    initial_holdup = _RawPhase(mass=5.0, mass_frac=(0.4, 0.6))
    unit = _UnitOperation(inlet)
    unit.__original_phase__ = initial_holdup
    sim = _sim_with_unit(unit)

    _, raw_cost, _ = sim.GetOPEX(
        1.0,
        include_holdups=False,
        kwargs_items={"raw_materials": {"basis": "mass"}},
    )

    assert "Initial_holdup" not in raw_cost.index.get_level_values(1)
    assert list(raw_cost.columns) == ["mass", "mass_A", "mass_B"]
    np.testing.assert_allclose(
        raw_cost.iloc[0][["mass", "mass_A", "mass_B"]],
        [20.0, 5.0, 15.0],
    )


def test_mole_basis_totals_use_canonical_singular_basis():
    inlet = _RawPhase(mole_flow=4.0, mole_frac=(0.25, 0.75))
    sim = _sim_with_unit(_UnitOperation(inlet))

    raw_materials = sim.GetRawMaterials(basis="mole")

    assert list(raw_materials.columns) == ["moles", "moles_A", "moles_B"]
    np.testing.assert_allclose(
        raw_materials.iloc[0][["moles", "moles_A", "moles_B"]],
        [40.0, 10.0, 30.0],
    )


def test_batch_raw_inlet_records_total_before_aggregation():
    inlet = _RawPhase(mass=6.0, mass_frac=(0.2, 0.8))
    unit = _UnitOperation(inlet)
    unit.oper_mode = "Batch"
    sim = _sim_with_unit(unit)

    raw_materials = sim.GetRawMaterials(basis="mass")

    assert list(raw_materials.columns) == ["mass", "mass_A", "mass_B"]
    np.testing.assert_allclose(
        raw_materials.iloc[0][["mass", "mass_A", "mass_B"]],
        [6.0, 1.2, 4.8],
    )


def test_mixed_phase_raw_inlet_is_decomposed_by_phase():
    liquid = _RawPhase(mass_flow=1.0, mass_frac=(0.8, 0.2))
    solid = _RawPhase(mass_flow=3.0, mass_frac=(0.1, 0.9))
    mixed = _MixedRawInlet([liquid, solid])
    sim = _sim_with_unit(_UnitOperation(mixed))

    raw_materials = sim.GetRawMaterials(basis="mass", totals=False)

    assert len(raw_materials) == 2
    np.testing.assert_allclose(
        np.sort(raw_materials["mass"].to_numpy()),
        [10.0, 30.0],
    )
