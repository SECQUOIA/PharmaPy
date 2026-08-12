"""Regression tests for DynamicCollector result bookkeeping (issue #75).

These tests drive the real ``DynamicCollector.solve_unit`` public path with real
inlet objects and real thermodynamic data, and stub only one expensive optional
collaborator per test so the checks stay in the core (non-Assimulo) lane.

The liquid-mixer test stubs the Assimulo integrator with a single explicit-Euler
step taken from the model's own right-hand side and initial state, so the
returned trajectory carries the production state layout rather than an assumed
one. The crystallizer test stubs the delegated ``SemibatchCryst`` sub-model and
asserts only how the collector dispatches plotting.
"""

import numpy as np
import pytest

import PharmaPy.Containers as containers
from PharmaPy.Containers import DynamicCollector
from PharmaPy.MixedPhases import SlurryStream
from PharmaPy.Streams import LiquidStream, SolidStream


pytestmark = pytest.mark.unit

INLET_MASS_FLOW = 20.0  # [kg/s]
# Asymmetric composition over the four species in pfr_test_pure_comp.json
# (A, B, C, solv) so a shifted or truncated split is visible element-wise.
INLET_MASS_FRAC = np.array([0.1, 0.2, 0.3, 0.4])  # [-]
INLET_TEMP = 320.0  # [K]
RUNTIME = 5.0  # [s]


class _RecordedProblem:
    """Stand-in for ``assimulo.problem.Explicit_Problem``.

    Records the right-hand side and initial state that ``solve_unit`` builds,
    without interpreting their layout.
    """

    def __init__(self, rhs, y0, t0=0.0):
        self.rhs = rhs
        self.y0 = np.asarray(y0, dtype=float)  # states_init, model layout
        self.t0 = t0  # [s]


class _EulerSolver:
    """Stand-in for ``assimulo.solvers.CVode`` taking one explicit-Euler step.

    The integrator is the only stubbed boundary; the derivative comes from the
    unit's own ``unit_model``, so the returned states keep whatever ordering
    the production code uses.
    """

    def __init__(self, problem):
        self.problem = problem
        self.verbosity = 0

    def simulate(self, final_time, ncp_list=None):
        """Return the two-point trajectory ``[y0, y0 + h * f(t0, y0)]``.

        Parameters
        ----------
        final_time : float
            End of the integration interval [s].
        ncp_list : array_like, optional
            Ignored; present for signature compatibility.

        Returns
        -------
        tuple of numpy.ndarray
            Times [s] with shape ``(2,)`` and states with shape ``(2, n)``.
        """
        y0 = self.problem.y0
        step = final_time - self.problem.t0  # [s]
        derivative = np.asarray(self.problem.rhs(self.problem.t0, y0),
                                dtype=float)
        states = np.vstack((y0, y0 + step * derivative))
        time = np.array([self.problem.t0, final_time])  # [s]

        return time, states


@pytest.fixture
def euler_backend(monkeypatch):
    """Replace the optional Assimulo constructors used by Containers.py."""
    monkeypatch.setattr(containers, "Explicit_Problem", _RecordedProblem)
    monkeypatch.setattr(containers, "CVode", _EulerSolver)


def test_liquid_mixer_requires_an_integration_end(data_path, euler_backend):
    """Reject a liquid solve with neither runtime nor requested times.

    Parameters
    ----------
    data_path : dict of pathlib.Path
        Repository test-data directories.
    euler_backend : None
        Fixture replacing the optional Assimulo integration boundary.
    """
    path = str(data_path["integration"] / "pfr_test_pure_comp.json")
    inlet = LiquidStream(path, temp=INLET_TEMP, mass_flow=INLET_MASS_FLOW,
                         mass_frac=INLET_MASS_FRAC)
    collector = DynamicCollector()
    collector.Inlet = inlet

    message = (r"DynamicCollector\.solve_unit requires 'runtime' \[s\] or "
               r"'time_grid' \[s\]; neither was supplied\.")
    with pytest.raises(ValueError, match=message):
        collector.solve_unit()


def test_liquid_mixer_result_labels_match_state_vector(data_path,
                                                       euler_backend):
    """``result`` must label holdup mass [kg] and composition [-] correctly.

    The holdup is created from the inlet stream, so its initial composition is
    the inlet composition and its temperature is the inlet temperature. With a
    constant inlet the only balance that moves is the total mass, whose exact
    solution over one step is ``mass_flow * step`` [kg].
    """
    path = str(data_path["integration"] / "pfr_test_pure_comp.json")

    inlet = LiquidStream(path, temp=INLET_TEMP, mass_flow=INLET_MASS_FLOW,
                         mass_frac=INLET_MASS_FRAC)

    collector = DynamicCollector()
    collector.Inlet = inlet

    time, _ = collector.solve_unit(runtime=RUNTIME)
    step = time[-1] - time[0]  # [s]

    result = collector.result

    # Composition is a mass fraction vector over the four species.
    assert result.mass_frac.shape == (2, len(INLET_MASS_FRAC))
    np.testing.assert_allclose(result.mass_frac[0], INLET_MASS_FRAC,
                               rtol=1e-10)
    np.testing.assert_allclose(result.mass_frac.sum(axis=1), np.ones(2),
                               rtol=1e-10)

    # Holdup mass accumulates the inlet mass flow: d(mass)/dt = mass_flow.
    assert result.mass.shape == (2,)
    assert result.mass[1] - result.mass[0] == pytest.approx(
        INLET_MASS_FLOW * step, rel=1e-10)

    # Temperature is unchanged because inlet and holdup enthalpies coincide.
    np.testing.assert_allclose(result.temp, np.full(2, INLET_TEMP),
                               rtol=1e-10)

    # `outputs` is the same labeled mapping handed downstream.
    np.testing.assert_allclose(collector.outputs["mass"], result.mass,
                               rtol=1e-10)
    np.testing.assert_allclose(collector.outputs["mass_frac"],
                               result.mass_frac, rtol=1e-10)


# --- Crystallizer-mode plotting dispatch -----------------------------------

SLURRY_VOL_FLOW = 1e-4  # [m**3/s]
SOLID_MASS_FRAC = np.array([0.0, 0.0, 1.0, 0.0])  # [-], solid is species C
TARGET_INDEX = 2  # index of species C in pfr_test_pure_comp.json
CRYST_SCALE = 1e-9  # [-], distribution scaling used by SemibatchCryst
CRYST_RUNTIME = 10.0  # [s]
NUM_CRYST_BINS = 15  # [-]


class _RecordingSemibatchCryst:
    """Stand-in for the delegated ``SemibatchCryst`` sub-model.

    Records the plotting delegation the collector is expected to perform and
    returns a state array shaped like a crystallizer solve, whose leading
    columns are crystal-size-distribution bins rather than liquid states.
    """

    #: Sentinels returned by :meth:`plot_profiles`.
    fig = object()
    axes = object()
    ax_right = object()

    def __init__(self, method=None, adiabatic=None, **kwargs):
        self.method = method
        self.adiabatic = adiabatic
        self.kwargs = kwargs

        self.elapsed_time = 0.0  # [s]
        self.states_di = {'distrib': {'dim': NUM_CRYST_BINS, 'type': 'diff'}}
        self.plot_calls = 0

        self.Outlet = _CrystOutlet()
        self.outputs = {}
        self.result = object()

    def solve_unit(self, runtime=None, time_grid=None, verbose=True):
        """Return a two-point CSD-leading trajectory without integrating.

        Parameters
        ----------
        runtime : float, optional
            Requested integration span [s].
        time_grid : array_like, optional
            Ignored.
        verbose : bool, optional
            Ignored.

        Returns
        -------
        tuple of numpy.ndarray
            Times [s] with shape ``(2,)`` and states with shape
            ``(2, NUM_CRYST_BINS + 5)``.
        """
        time = np.array([self.elapsed_time,
                         self.elapsed_time + runtime])  # [s]
        num_columns = NUM_CRYST_BINS + len(SOLID_MASS_FRAC) + 1
        states = np.tile(np.arange(num_columns, dtype=float), (2, 1))

        return time, states

    def plot_profiles(self, fig_size=None, time_div=1, **kwargs):
        """Record the delegated call and return sentinel figure handles."""
        self.plot_calls += 1

        return self.fig, self.axes, self.ax_right


class _CrystOutlet:
    """Minimal outlet exposing the volume the collector reads back."""

    vol = 1e-3  # [m**3]


@pytest.fixture
def slurry_inlet(data_path):
    """Real crystallizing inlet built from the repository thermo data.

    A ``DynamicCollector`` only ever receives a crystallizing inlet from an
    upstream unit, so ``y_inlet``/``y_upstream`` are populated here the way
    ``PharmaPy.Connections.Connection`` populates them; without them the inlet
    concentration would fall back to the zero-filled missing-state default.
    """
    path = str(data_path["integration"] / "pfr_test_pure_comp.json")

    x_distrib = np.geomspace(1.0, 1500.0, num=NUM_CRYST_BINS)  # [um]
    # Narrow seed distribution centred on the grid [#/m**3/um].
    distrib = np.exp(-((np.log(x_distrib) - np.log(100.0)) / 0.5)**2) * 1e10

    liquid = LiquidStream(path, temp=310.0, mass_frac=[0.1, 0.2, 0.3, 0.4],
                          vol_flow=SLURRY_VOL_FLOW)
    solid = SolidStream(path, temp=310.0, x_distrib=x_distrib, distrib=distrib,
                        mass_frac=SOLID_MASS_FRAC)

    inlet = SlurryStream(vol_flow=SLURRY_VOL_FLOW, x_distrib=x_distrib,
                         distrib=distrib)
    inlet.Phases = [liquid, solid]

    inlet.y_inlet = {'mass_conc': liquid.mass_conc,  # [kg/m**3]
                     'vol_flow': SLURRY_VOL_FLOW,  # [m**3/s]
                     'temp': inlet.temp,  # [K]
                     'distrib': distrib}  # [#/m**3/um]
    inlet.y_upstream = inlet.y_inlet
    inlet.time_upstream = None

    return inlet


def test_inlet_assignment_sets_collector_model_mode(data_path, slurry_inlet):
    """Derive and refresh collector mode whenever an inlet is assigned.

    Parameters
    ----------
    data_path : dict of pathlib.Path
        Repository test-data directories.
    slurry_inlet : SlurryStream
        Real crystallizing inlet with composition [-], flow [m**3/s],
        temperature [K], and distribution [#/m**3/um] states.
    """
    collector = DynamicCollector()

    assert collector.is_cryst is False

    collector.Inlet = slurry_inlet

    assert collector.is_cryst is True

    path = str(data_path["integration"] / "pfr_test_pure_comp.json")
    liquid_inlet = LiquidStream(path, temp=INLET_TEMP,
                                mass_flow=INLET_MASS_FLOW,
                                mass_frac=INLET_MASS_FRAC)
    collector.Inlet = liquid_inlet

    assert collector.is_cryst is False


def test_crystallizer_collector_delegates_plotting(slurry_inlet, monkeypatch):
    """A crystallizing collector must plot through its crystallizer sub-model.

    ``plot_local`` reads liquid-shaped slices of the state vector, which for a
    crystallizer solve are crystal-size-distribution bins. After a crystallizer
    solve, ``plot_profiles`` must therefore delegate to the sub-model's own
    plotter instead.
    """
    monkeypatch.setattr(containers, "SemibatchCryst", _RecordingSemibatchCryst)

    collector = DynamicCollector()
    collector.Inlet = slurry_inlet
    collector.KinCryst = object()
    collector.kwargs_cryst = {'target_ind': TARGET_INDEX, 'target_comp': 'C',
                              'scale': CRYST_SCALE}

    collector.solve_unit(runtime=CRYST_RUNTIME)

    assert collector.model_type == 'crystallizer'

    fig, axes = collector.plot_profiles()

    assert collector.CrystInst.plot_calls == 1
    assert fig is _RecordingSemibatchCryst.fig
    assert axes is _RecordingSemibatchCryst.axes


def test_crystallizer_results_skip_liquid_profile_slicing(slurry_inlet,
                                                          monkeypatch):
    """Crystallizer retrieval must not overwrite liquid-profile attributes.

    Parameters
    ----------
    slurry_inlet : SlurryStream
        Real crystallizing inlet with composition [-], flow [m**3/s],
        temperature [K], and distribution [#/m**3/um] states.
    monkeypatch : pytest.MonkeyPatch
        Fixture used to replace only the solver-backed crystallizer boundary.
    """
    monkeypatch.setattr(containers, "SemibatchCryst", _RecordingSemibatchCryst)

    collector = DynamicCollector()
    collector.Inlet = slurry_inlet
    collector.KinCryst = object()
    collector.kwargs_cryst = {'target_ind': TARGET_INDEX, 'target_comp': 'C',
                              'scale': CRYST_SCALE}
    liquid_profile_sentinel = object()
    collector.wConcProf = liquid_profile_sentinel
    collector.massProf = liquid_profile_sentinel
    collector.tempProf = liquid_profile_sentinel

    collector.solve_unit(runtime=CRYST_RUNTIME)

    assert collector.wConcProf is liquid_profile_sentinel
    assert collector.massProf is liquid_profile_sentinel
    assert collector.tempProf is liquid_profile_sentinel
