"""Aggregation contracts for MixedPhase and MixedStream.

MixedPhase has no explicit accessors for most quantities. Attribute lookups
fall through ``__getattr__``, which sums the members when the name is listed
in ``EXTENSIVE_PROPERTIES`` and otherwise returns a mass-weighted average.
Getting that classification wrong is silent, so these tests pin both branches.

The fixture is deliberately asymmetric: a 3 kg liquid against a 1 kg solid, so
that a sum (4 kg) and a mass-weighted average (2.5 kg) are far apart and
cannot be confused, and phase temperatures of 310 K and 290 K give a weighted
mean of 305 K that no unweighted mean would reproduce.
"""

import os

import numpy as np
import pytest

from PharmaPy.IntegratorBackends import AssimuloBackend
from PharmaPy.MixedPhases import MixedPhase, MixedStream
from PharmaPy.Phases import LiquidPhase, SolidPhase
from PharmaPy.ProcessControl_Refactor import DefaultContinuousVesselVolume
from PharmaPy.Reactors_refactor import ContinuousReactor
from PharmaPy.Streams import LiquidStream, SolidStream

pytestmark = pytest.mark.unit

DATA_PATH = os.path.join(
    os.path.dirname(__file__), "Flowsheet", "data", "compound_database.json"
)

LIQUID_MASS = 3.0  # [kg]
SOLID_MASS = 1.0  # [kg]
LIQUID_MASS_FRAC = [0.4, 0.6, 0.0, 0.0, 0.0]  # [-]
SOLID_MASS_FRAC = [0.0, 0.0, 1.0, 0.0, 0.0]  # [-]
LIQUID_TEMP = 310.0  # [K]
SOLID_TEMP = 290.0  # [K]

LIQUID_MASS_FLOW = 0.3  # [kg/s]
SOLID_MASS_FLOW = 0.1  # [kg/s]

VESSEL_DIAMETER = 0.01  # [m]

# Aggregates are compared against sums of the same floating point members, so
# only accumulated roundoff is at stake.
AGGREGATE_RTOL = 1e-12  # [-]


def _liquid(mass=LIQUID_MASS, temp=LIQUID_TEMP):
    """Build the liquid member of the fixture.

    Parameters
    ----------
    mass : float, optional
        Liquid holdup [kg].
    temp : float, optional
        Liquid temperature [K].

    Returns
    -------
    LiquidPhase
        Liquid phase carrying the standard A/B composition.
    """
    return LiquidPhase(
        DATA_PATH, mass=mass, mass_frac=LIQUID_MASS_FRAC, temp=temp
    )


def _solid(mass=SOLID_MASS, temp=SOLID_TEMP):
    """Build the solid member of the fixture.

    Parameters
    ----------
    mass : float, optional
        Solid holdup [kg].
    temp : float, optional
        Solid temperature [K].

    Returns
    -------
    SolidPhase
        Solid phase of pure species C.
    """
    return SolidPhase(
        DATA_PATH, mass=mass, mass_frac=SOLID_MASS_FRAC, temp=temp
    )


def _streams():
    """Build the two-phase feed members.

    Returns
    -------
    tuple of (LiquidStream, SolidStream)
        Liquid and solid streams at ``LIQUID_MASS_FLOW`` and
        ``SOLID_MASS_FLOW`` [kg/s].
    """
    liquid = LiquidStream(
        DATA_PATH, mass_flow=LIQUID_MASS_FLOW, mass_frac=LIQUID_MASS_FRAC
    )
    solid = SolidStream(
        DATA_PATH, mass_flow=SOLID_MASS_FLOW, mass_frac=SOLID_MASS_FRAC
    )
    return liquid, solid


def test_mixed_phase_sums_extensive_amounts():
    """Holdup quantities add across the constituent phases."""
    liquid, solid = _liquid(), _solid()
    mixed = MixedPhase([liquid, solid])

    np.testing.assert_allclose(
        mixed.mass, LIQUID_MASS + SOLID_MASS, rtol=AGGREGATE_RTOL
    )
    np.testing.assert_allclose(
        mixed.vol, liquid.vol + solid.vol, rtol=AGGREGATE_RTOL
    )
    np.testing.assert_allclose(
        mixed.moles, liquid.moles + solid.moles, rtol=AGGREGATE_RTOL
    )


def test_mixed_phase_mass_weights_intensive_properties():
    """Temperature is a mass-weighted mean, not a sum.

    Guards the opposite failure to the extensive case: classifying an
    intensive property as extensive would report a physically impossible
    600 K for two phases near ambient.
    """
    mixed = MixedPhase([_liquid(), _solid()])

    total_mass = LIQUID_MASS + SOLID_MASS  # [kg]
    expected = (
        LIQUID_MASS * LIQUID_TEMP + SOLID_MASS * SOLID_TEMP
    ) / total_mass  # [K], equals 305.0
    np.testing.assert_allclose(mixed.temp, expected, rtol=AGGREGATE_RTOL)


def test_mixed_stream_sums_flow_rates():
    """Flow rates add across the phases of a multiphase stream.

    A mass-weighted average would report 2.5 kg/s for a 3 kg/s liquid and a
    1 kg/s solid, understating the true 4 kg/s throughput.
    """
    liquid, solid = _streams()
    stream = MixedStream([liquid, solid])

    np.testing.assert_allclose(
        stream.mass_flow,
        LIQUID_MASS_FLOW + SOLID_MASS_FLOW,
        rtol=AGGREGATE_RTOL,
    )
    np.testing.assert_allclose(
        stream.vol_flow,
        liquid.vol_flow + solid.vol_flow,
        rtol=AGGREGATE_RTOL,
    )
    np.testing.assert_allclose(
        stream.mole_flow,
        liquid.mole_flow + solid.mole_flow,
        rtol=AGGREGATE_RTOL,
    )


def test_mixed_stream_sums_amount_aliases():
    """A stream still adds the amount names that alias its flow rates.

    Stream objects keep ``mass`` as an alias of ``mass_flow``, so the stream
    classification must extend the phase classification rather than replace
    it. Replacing it fixes the flow names while silently turning ``mass``
    into an average.
    """
    liquid, solid = _streams()
    stream = MixedStream([liquid, solid])

    np.testing.assert_allclose(
        stream.mass, liquid.mass + solid.mass, rtol=AGGREGATE_RTOL
    )


def test_stream_built_by_conversion_sums_flow_rates():
    """``to_stream`` produces a stream with stream aggregation rules.

    The conversion copies the instance and reassigns ``__class__`` without
    running ``MixedStream.__init__``, so a per-instance classification set on
    the phase would survive into the stream.
    """
    liquid, solid = _liquid(), _solid()
    stream = MixedPhase([liquid, solid]).to_stream()

    assert isinstance(stream, MixedStream)
    np.testing.assert_allclose(
        stream.mass_flow, LIQUID_MASS + SOLID_MASS, rtol=AGGREGATE_RTOL
    )


def test_continuous_vessel_withdraws_the_total_fed_volume():
    """The volume controller sees the whole two-phase feed.

    DefaultContinuousVesselVolume reads the aggregate ``vol_flow`` of each
    inlet stream to hold the vessel level. Averaging instead of summing makes
    the controller withdraw less than is fed, so the vessel silently fills.
    """
    # ContinuousReactor takes its controller from a mutable default argument,
    # so every vessel built without an explicit one shares a single instance
    # and inherits whichever target volume was latched first. Pass a fresh
    # controller until that default is fixed.
    vessel = ContinuousReactor(
        integrator=AssimuloBackend(),
        h_conv=0,
        diam=VESSEL_DIAMETER,
        isothermal=True,
        controller=DefaultContinuousVesselVolume(),
    )
    vessel.Phases = [_liquid(temp=298.15), _solid(temp=298.15)]

    liquid, solid = _streams()
    vessel.Inlet = MixedStream([liquid, solid])

    states = vessel.create_solver_init_states()
    vessel.unit_model(0.0, states, mat_bce=True)

    expected_vol_flow = liquid.vol_flow + solid.vol_flow  # [m**3/s]
    posted = [
        value
        for key, value in vessel.controller.operating_conditions.items()
        if key.name == "vol_flow" and key.port == "outlet"
    ]

    assert len(posted) == 1
    np.testing.assert_allclose(
        posted[0], expected_vol_flow, rtol=AGGREGATE_RTOL
    )
