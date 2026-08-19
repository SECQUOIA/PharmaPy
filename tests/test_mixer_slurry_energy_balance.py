"""Regression tests for the Mixer adiabatic energy balance with solids.

Issue #25: ``Mixer.energy_balance`` weights every inlet enthalpy by the
liquid-only mass ``u_inputs['mass_liq']`` and closes the balance against
``sum(mass_liq) * h_out``, even though a solids-bearing inlet reports a
mixed-phase enthalpy. The solid mass therefore carries no enthalpy on either
side of the balance and the solved adiabatic mixing temperature is wrong
whenever a solid phase is present.

The same lines also select the wrong enthalpy basis: the inlet call is shaped
for ``Cake.getEnthalpy(temp, mass_frac=..., distrib=...)`` and raises
``TypeError`` for a ``Slurry``/``SlurryStream`` inlet, while the outlet call
``self.Outlet.getEnthalpy(temp)`` takes ``Slurry.getEnthalpy``'s default
``volumetric=True`` and returns [J/m**3] of suspension, which is then
multiplied by a mass in [kg]. A correct balance must use the mass-basis
mixed-phase enthalpy ``getEnthalpy(temp, volumetric=False)`` [J/kg] on both
sides and weight it by the total slurry mass ``mass_liq + mass_solid`` [kg].

Fixtures use the shared ``tests/Flowsheet/data/compound_database.json`` thermo
file and real ``LiquidPhase``/``SolidPhase``/``Slurry`` collaborators. Expected
temperatures are derived independently of the production expression, by closing
the total-enthalpy balance over the individual phase enthalpies with ``brentq``.

Both tests call ``Mixer.energy_balance`` with the ``u_inputs`` mapping produced
by the real ``Mixer.get_inputs_solids``, and assign ``Mixer.Outlet`` the way
``Mixer.balances_solids`` does. They cannot yet run through the public
``Mixer.solve_unit`` entry point because two separate, unfixed defects abort
that path first, both recorded in the triage comment on issue #25 and not
covered by #38 or #88:

* ``Containers.py:527`` constructs ``LiquidPhase(path)`` with no composition,
  raising ``ValueError: No measure of composition was provided`` for every
  solids-bearing inlet set.
* the outlet ``SolidPhase(path, mass=..., distrib=...)`` built at
  ``Containers.py:487-497`` yields ``nan`` moments, so assigning
  ``Slurry.Phases`` raises ``RuntimeError`` from ``newton``.

The provisional contract asserted here is therefore the energy balance alone.
Once those blockers are fixed, these tests should drive ``Mixer.solve_unit``
end to end and assert the outlet temperature it reports.
"""

import numpy as np
import pytest
from scipy.optimize import brentq

from PharmaPy.Containers import Mixer
from PharmaPy.MixedPhases import Slurry
from PharmaPy.Phases import LiquidPhase, SolidPhase


pytestmark = pytest.mark.unit

# Liquid composition [-] over the database species (A, B, C, D, solvent).
# Identical in both inlets so that the mixed liquid composition, and hence the
# liquid density, is unchanged by mixing and volumes are exactly additive.
LIQUID_MASS_FRAC = [0.1, 0.1, 0.1, 0.1, 0.6]

# Solid composition [-]: pure component A.
SOLID_MASS_FRAC = [1.0, 0.0, 0.0, 0.0, 0.0]

# Volumetric shape factor of the crystals [-].
KV = 1.0

# Crystal size grid [um] and a volume-specific number density [#/m**3/um]
# whose third moment is 0.2 m**3/m**3 by the trapezoidal rule on the 100 um
# grid: 100 * (1e9 * 100**3 + 1.25e8 * 200**3) * 1e-18 = 0.2.
X_DISTRIB = np.array([0.0, 100.0, 200.0, 300.0])  # [um]
DISTRIB = np.array([0.0, 1.0e9, 1.25e8, 0.0])  # [#/m**3/um]

# Mixer geometry-free operating point.
MASS_LIQUID_INLET = 1.0  # [kg], solids-free inlet
VOL_SLURRY_INLET = 1.0e-3  # [m**3], slurry inlet

# Bracket for the independent outlet-temperature solve [K]. Wide enough to
# contain any physically admissible adiabatic mixing temperature for the
# inlet temperatures used below.
TEMP_BRACKET = (250.0, 450.0)  # [K]


@pytest.fixture(scope="module")
def thermo_path(data_path):
    """Path of the shared pure-component thermodynamic database."""
    return str(data_path["flowsheet"] / "compound_database.json")


def _make_slurry_inlet(thermo_path, temp):
    """Build a slurry inlet at a uniform temperature.

    Parameters
    ----------
    thermo_path : str
        Path to the pure-component thermodynamic database.
    temp : float
        Temperature of both slurry phases [K].

    Returns
    -------
    PharmaPy.MixedPhases.Slurry
        Slurry of volume ``VOL_SLURRY_INLET`` [m**3] whose liquid and solid
        masses [kg] follow from the volume-specific distribution.
    """
    liquid = LiquidPhase(thermo_path, mass_frac=LIQUID_MASS_FRAC, temp=temp)
    solid = SolidPhase(thermo_path, mass_frac=SOLID_MASS_FRAC, kv=KV,
                       temp=temp)

    slurry = Slurry(vol=VOL_SLURRY_INLET, x_distrib=X_DISTRIB,
                    distrib=DISTRIB)
    slurry.Phases = (liquid, solid)

    return slurry


def _build_mixer(thermo_path, temp_liquid, temp_slurry):
    """Assemble a Mixer fed one solids-free inlet and one slurry inlet.

    Parameters
    ----------
    thermo_path : str
        Path to the pure-component thermodynamic database.
    temp_liquid : float
        Temperature of the solids-free liquid inlet [K].
    temp_slurry : float
        Temperature of the slurry inlet [K].

    Returns
    -------
    mixer : PharmaPy.Containers.Mixer
        Mixer with both inlets attached and ``Outlet`` set to the mixed
        slurry, as ``Mixer.balances_solids`` does before solving the energy
        balance.
    u_inputs : dict
        Inlet mapping returned by the real ``Mixer.get_inputs_solids``.
    masses : dict
        Total liquid and solid mass of the mixture [kg], keyed ``'liquid'``
        and ``'solid'``.
    """
    liquid_inlet = LiquidPhase(thermo_path, mass=MASS_LIQUID_INLET,
                               mass_frac=LIQUID_MASS_FRAC, temp=temp_liquid)
    slurry_inlet = _make_slurry_inlet(thermo_path, temp_slurry)

    mixer = Mixer()
    # The solids-free inlet is listed first because ``Mixer.Inlets`` reads
    # ``name_species`` off ``Inlets[0]``, which ``Slurry`` does not define.
    mixer.Inlets = [liquid_inlet, slurry_inlet]

    u_inputs, _ = mixer.get_inputs_solids()

    total_liquid = float(u_inputs['mass_liq'].sum())  # [kg]
    total_solid = float(u_inputs['mass_solid'].sum())  # [kg]

    # Additive volumes: both inlets share the liquid composition, so the
    # mixed liquid density equals each inlet's liquid density.
    dens_liquid = slurry_inlet.Liquid_1.getDensity(
        mass_frac=np.array(LIQUID_MASS_FRAC))  # [kg/m**3]
    vol_liquid_inlet = MASS_LIQUID_INLET / dens_liquid  # [m**3]
    vol_total = vol_liquid_inlet + VOL_SLURRY_INLET  # [m**3]

    # Volume-specific outlet distribution [#/m**3/um]: the solids-free inlet
    # contributes no crystals, so the slurry population is simply diluted.
    distrib_out = DISTRIB * VOL_SLURRY_INLET / vol_total

    outlet = Slurry(vol=vol_total, x_distrib=X_DISTRIB, distrib=distrib_out)
    outlet.Phases = (
        LiquidPhase(thermo_path, mass_frac=LIQUID_MASS_FRAC),
        SolidPhase(thermo_path, mass_frac=SOLID_MASS_FRAC, kv=KV),
    )

    # Fixture self-check: the constructed outlet must hold the same inventory
    # the material balance produces, otherwise the energy comparison below is
    # meaningless.
    assert outlet.Liquid_1.mass == pytest.approx(total_liquid, rel=1e-10)
    assert outlet.Solid_1.mass == pytest.approx(total_solid, rel=1e-10)

    mixer.Outlet = outlet

    return mixer, u_inputs, {'liquid': total_liquid, 'solid': total_solid}


def _expected_outlet_temp(mixer, u_inputs, masses):
    """Solve the adiabatic outlet temperature from the phase enthalpies.

    The balance is closed directly over the individual liquid and solid phase
    enthalpies [J/kg], independently of ``Mixer.energy_balance`` and of
    ``Slurry.getEnthalpy``::

        sum_i (m_liq_i h_liq(T_i) + m_sol_i h_sol(T_i))
            = M_liq h_liq(T_out) + M_sol h_sol(T_out)

    Parameters
    ----------
    mixer : PharmaPy.Containers.Mixer
        Mixer whose inlets supply the phase collaborators.
    u_inputs : dict
        Inlet mapping from ``Mixer.get_inputs_solids``.
    masses : dict
        Total liquid and solid mass of the mixture [kg].

    Returns
    -------
    float
        Adiabatic outlet temperature [K].
    """
    liquid_probe = mixer.Outlet.Liquid_1
    solid_probe = mixer.Outlet.Solid_1

    enthalpy_in = 0.0  # [J]
    for ind, temp in enumerate(u_inputs['temp']):
        enthalpy_in += u_inputs['mass_liq'][ind] * liquid_probe.getEnthalpy(
            temp=temp, mass_frac=u_inputs['mass_frac'][ind], basis='mass')
        enthalpy_in += u_inputs['mass_solid'][ind] * solid_probe.getEnthalpy(
            temp=temp, basis='mass')

    def residual(temp):
        enthalpy_out = (
            masses['liquid'] * liquid_probe.getEnthalpy(
                temp=temp, mass_frac=np.array(LIQUID_MASS_FRAC), basis='mass')
            + masses['solid'] * solid_probe.getEnthalpy(temp=temp,
                                                        basis='mass')
        )  # [J]

        return enthalpy_in - enthalpy_out

    return brentq(residual, *TEMP_BRACKET, xtol=1e-12, rtol=1e-14)


def test_mixer_energy_balance_preserves_isothermal_inlet_temperature(
        thermo_path):
    """Mixing streams already at a common temperature must not shift it.

    An adiabatic mixer fed a solids-free liquid and a slurry, both at
    ``temp_common``, must return ``temp_common``: no enthalpy is exchanged.
    Weighting the mixed-phase inlet enthalpy by liquid mass while the outlet
    carries the solid mass as well breaks this invariant.
    """
    temp_common = 300.0  # [K]

    mixer, u_inputs, _ = _build_mixer(thermo_path, temp_common, temp_common)

    temp_out = mixer.energy_balance(u_inputs)  # [K]

    assert temp_out == pytest.approx(temp_common, abs=1e-8)


def test_mixer_energy_balance_weights_inlets_by_total_slurry_mass(
        thermo_path):
    """The solved temperature must close a total-slurry-mass enthalpy balance.

    A hot solids-free inlet and a cold slurry inlet have different solid-to-
    liquid mass ratios, so a liquid-only weighting cannot be absorbed into a
    common factor and the solved temperature differs from the correct one.
    """
    temp_liquid = 350.0  # [K]
    temp_slurry = 300.0  # [K]

    mixer, u_inputs, masses = _build_mixer(thermo_path, temp_liquid,
                                           temp_slurry)

    expected_temp = _expected_outlet_temp(mixer, u_inputs, masses)  # [K]

    # Sanity: an adiabatic mix must land strictly between the inlet
    # temperatures, so the reference value cannot be trivially satisfied.
    assert temp_slurry < expected_temp < temp_liquid

    temp_out = mixer.energy_balance(u_inputs)  # [K]

    assert temp_out == pytest.approx(expected_temp, rel=1e-9)
