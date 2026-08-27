"""Particle-size unit regressions for deliquoring and drying setup.

``SolidPhase.x_distrib`` stores crystal sizes in micrometers, while the
capillary and threshold-pressure correlations consume particle diameters in
meters. These tests stop both units at the solver boundary so the assertions
exercise the real setup calculations without running an Assimulo transient.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from PharmaPy.SolidLiquidSep import DeliquoringStep
from PharmaPy import SolidLiquidSep as solid_liquid_sep
from PharmaPy import Drying_Model as drying_model


pytestmark = pytest.mark.unit


SIZE_GRID_UM = np.array([50.0, 100.0, 150.0])  # [um]
CSD_NUMBER = np.array([1.0, 2.0, 1.0])  # [#/m**3/um]


class StopAfterSetup(Exception):
    """Stop a setup-only solve before time integration begins."""


class DummyExplicitProblem:
    """Minimal ``Explicit_Problem`` replacement for setup-only tests."""

    def __init__(self, rhs, y0, t0, **kwargs):
        """Store solver problem inputs without constructing Assimulo objects.

        Parameters
        ----------
        rhs : callable
            Right-hand-side function that would be integrated [-].
        y0 : ndarray
            Initial flattened state passed to the solver [-].
        t0 : float
            Initial integration time [-].
        **kwargs : dict
            Additional problem options, such as event switches [-].
        """
        self.rhs = rhs
        self.y0 = y0
        self.t0 = t0
        self.kwargs = kwargs
        self.name = None


class StopSolver:
    """Minimal ``CVode`` replacement that stops at ``simulate``."""

    def __init__(self, problem):
        """Store the explicit problem supplied by ``solve_unit``.

        Parameters
        ----------
        problem : DummyExplicitProblem
            Problem object constructed after setup calculations complete [-].
        """
        self.problem = problem
        self.linear_solver = None
        self.verbosity = None

    def simulate(self, *args, **kwargs):
        """Stop before numerical integration.

        Parameters
        ----------
        *args : tuple
            Positional solver arguments that would include final time [-].
        **kwargs : dict
            Keyword solver arguments [-].

        Raises
        ------
        StopAfterSetup
            Always raised to keep the test focused on setup calculations.
        """
        raise StopAfterSetup


def _zeroth_moment() -> float:
    """Return the number-based CSD zeroth moment.

    Returns
    -------
    float
        Integral of ``CSD_NUMBER`` [#/m**3/um] over ``SIZE_GRID_UM`` [um],
        giving the total crystal count per suspension volume [#/m**3].
    """
    return solid_liquid_sep.trapezoidal_rule(SIZE_GRID_UM, CSD_NUMBER)


def test_deliquoring_setup_converts_micrometer_grid_to_meter_diameters(
    monkeypatch,
):
    """Deliquoring capillary quantities use meter diameters from ``x_distrib``."""
    captured = {}  # [-]
    size_grid_m = SIZE_GRID_UM * 1e-6  # [m]
    porosity = 0.45  # [-]
    surface_tension = np.array([0.065, 0.067])  # [N/m]
    surface_tension_mean = np.mean(surface_tension)  # [N/m]

    def capture_sat_inf(x_vec, csd, deltaP, porosity, height, mu_zero, props):
        """Capture the diameter grid supplied to ``get_sat_inf``.

        Parameters
        ----------
        x_vec : ndarray
            Particle diameter grid supplied to the capillary correlation [m].
        csd : ndarray
            Crystal-size distribution on the stored grid [#/m**3/um].
        deltaP : float
            Pressure drop supplied to the saturation correlation [Pa].
        porosity : float
            Cake porosity [-].
        height : float
            Cake height [m].
        mu_zero : float
            Zeroth distribution moment [#/m**3].
        props : tuple
            Surface tension [N/m] and liquid density [kg/m**3].

        Returns
        -------
        float
            Synthetic irreducible saturation [-].
        """
        captured["x_vec"] = x_vec.copy()  # [m]
        return 0.25  # [-]

    monkeypatch.setattr(solid_liquid_sep, "get_sat_inf", capture_sat_inf)
    monkeypatch.setattr(solid_liquid_sep, "Explicit_Problem", DummyExplicitProblem)
    monkeypatch.setattr(solid_liquid_sep, "CVode", StopSolver)

    unit = DeliquoringStep(num_nodes=2)
    unit.cake_height = 0.02  # [m]
    unit.z_centers = np.array([0.25, 0.75])  # [-]
    unit.delta_z = np.array([0.5, 0.5])  # [-]
    unit.Solid_1 = SimpleNamespace(
        distrib=CSD_NUMBER,
        x_distrib=SIZE_GRID_UM,
        moments=np.array([_zeroth_moment()]),  # [#/m**3]
        getPorosity=lambda: porosity,
        getDensity=lambda: 1500.0,  # [kg/m**3]
    )
    unit.Liquid_1 = SimpleNamespace(
        num_species=2,
        name_species=["solute_a", "solute_b"],
        getDensity=lambda: np.array([950.0, 970.0]),  # [kg/m**3]
        getViscosity=lambda: np.array([1.0e-3, 1.2e-3]),  # [Pa*s]
        getSurfTension=lambda: surface_tension,
        getDensityPure=lambda: np.array([1000.0, 1100.0]),  # [kg/m**3]
    )
    unit.CakePhase = SimpleNamespace(
        alpha=2.0e10,  # [m/kg]
        saturation=np.array([0.8]),  # [-]
        Liquid_1=SimpleNamespace(mass_conc=np.array([1.0, 2.0])),  # [kg/m**3]
        z_external=np.array([0.0, unit.cake_height]),  # [m]
    )

    with pytest.raises(StopAfterSetup):
        unit.solve_unit(deltaP=5.0e4, runtime=10.0)  # [Pa], [s]

    np.testing.assert_allclose(captured["x_vec"], size_grid_m)

    threshold_by_size = (
        4.6 * (1 - porosity) * surface_tension_mean / porosity / size_grid_m
    )  # [Pa]
    expected_threshold = (
        solid_liquid_sep.trapezoidal_rule(
            SIZE_GRID_UM,
            threshold_by_size * CSD_NUMBER,
        )
        / unit.Solid_1.moments[0]
    )  # [Pa]

    assert unit.p_thresh == pytest.approx(expected_threshold)


def test_drying_setup_converts_micrometer_grid_before_irreducible_saturation(
    monkeypatch,
):
    """Drying passes meter diameters to the shared saturation correlation."""
    captured = {}  # [-]

    def capture_sat_inf(x_vec, csd, deltaP, porosity, height, mu_zero, props):
        """Capture the diameter grid supplied by ``Drying.solve_unit``.

        Parameters
        ----------
        x_vec : ndarray
            Particle diameter grid supplied to the capillary correlation [m].
        csd : ndarray
            Crystal-size distribution on the stored grid [#/m**3/um].
        deltaP : float
            Pressure drop after medium resistance is removed [Pa].
        porosity : float
            Cake porosity [-].
        height : float
            Cake height [m].
        mu_zero : float
            Zeroth distribution moment [#/m**3].
        props : tuple
            Mean surface tension [N/m] and liquid density [kg/m**3].

        Returns
        -------
        float
            Synthetic irreducible saturation [-].
        """
        captured["x_vec"] = x_vec.copy()  # [m]
        return 0.25  # [-]

    monkeypatch.setattr(drying_model, "get_sat_inf", capture_sat_inf)

    dryer = drying_model.Drying(number_nodes=2, supercrit_names=["nitrogen"])
    dryer.names_states_in = ["temp", "mass_frac"]
    dryer.idx_supercrit = np.array([1])  # component indices [-]
    dryer.cake_height = 0.02  # [m]
    dryer.elapsed_time = 0.0  # [s]
    dryer.Liquid_1 = SimpleNamespace(
        num_species=3,
        getDensity=lambda temp, mass_frac, basis: np.array([930.0, 940.0]),  # [kg/m**3]
        getSurfTension=lambda temp, mass_frac: np.array([0.068, 0.070]),  # [N/m]
    )
    solid = SimpleNamespace(
        temp=302.0,  # [K]
        x_distrib=SIZE_GRID_UM,
        distrib=CSD_NUMBER,
        moments=np.array([_zeroth_moment()]),  # [#/m**3]
        getDensity=lambda: 1500.0,  # [kg/m**3]
        getCp=lambda: 700.0,  # [J/kg/K]
        getMoments=lambda mom_num: np.array([1.0, 1.0e-4, 1.0e-8, 1.0e-12, 0.0]),
    )
    dryer.Solid_1 = solid
    dryer.Vapor_1 = SimpleNamespace(
        mass_frac=np.array([0.01, 0.98, 0.01]),  # [-]
        temp=300.0,  # [K]
    )
    dryer.CakePhase = SimpleNamespace(
        Liquid_1=SimpleNamespace(mass_frac=np.array([0.20, 0.10, 0.70])),  # [-]
        Solid_1=solid,
        saturation=np.array([0.55]),  # [-]
        z_external=np.array([0.0, dryer.cake_height]),  # [m]
        alpha=2.0e10,  # [m/kg]
        porosity=0.45,  # [-]
    )

    def stop_at_rhs(time, states, sw=None):
        """Stop the drying setup at the RHS evaluation boundary.

        Parameters
        ----------
        time : float
            Initial solve time [s].
        states : ndarray
            Flattened initial drying state [-] and [K].
        sw : list, optional
            Assimulo event switches [-].

        Raises
        ------
        StopAfterSetup
            Always raised to keep the test focused on setup calculations.
        """
        raise StopAfterSetup

    monkeypatch.setattr(dryer, "unit_model", stop_at_rhs)

    with pytest.raises(StopAfterSetup):
        dryer.solve_unit(deltaP=5.0e4, runtime=10.0)  # [Pa], [s]

    np.testing.assert_allclose(captured["x_vec"], SIZE_GRID_UM * 1e-6)
