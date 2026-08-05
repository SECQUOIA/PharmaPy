"""Regression tests for the shape factor used by ``Cake.get_alpha``.

Issue #160 records that the scalar volumetric shape factor cancels from the
normalized crystal-volume weights. The hydraulic resistance ``alpha`` should
therefore be independent of any positive scalar ``kv``, while the method must
still obtain that factor from its attached solid phase rather than a literal.
"""

from unittest.mock import MagicMock, PropertyMock, patch

import numpy as np
import pytest

from PharmaPy.MixedPhases import Cake
from PharmaPy.Phases import LiquidPhase, SolidPhase


pytestmark = pytest.mark.unit

# Historical dimensionless literal replaced by the phase-owned value.
LEGACY_SHAPE_FACTOR = 0.524  # [-]

# Additional positive probes: 0.5 reproduces the #158 cake case, while 0.8 is
# a distinct scalar selected to exercise cancellation rather than calibration.
SHAPE_FACTOR_PROBES = (LEGACY_SHAPE_FACTOR, 0.5, 0.8)  # [-]

# Alpha recorded for the #158 fixture on unchanged master with kv = 0.524.
EXPECTED_LEGACY_ALPHA = 3.858303591823321e7  # [m/kg]

# Liquid composition for the shared A/B/C/D/solvent database [mass fraction].
LIQUID_MASS_FRAC = [0.1, 0.1, 0.1, 0.1, 0.6]  # [-]

# Pure-component-A solid composition [mass fraction].
SOLID_MASS_FRAC = [1.0, 0.0, 0.0, 0.0, 0.0]  # [-]


@pytest.fixture(scope="module")
def thermo_path(data_path):
    """Return the shared pure-component thermodynamic database path.

    Parameters
    ----------
    data_path : dict[str, pathlib.Path]
        Repository test-data directories.

    Returns
    -------
    str
        Path to the flowsheet compound database.
    """
    return str(data_path["flowsheet"] / "compound_database.json")


def _make_cake(thermo_path: str, shape_factor: float) -> Cake:
    """Build the #158 cake fixture with a selected crystal shape factor.

    Parameters
    ----------
    thermo_path : str
        Path to the pure-component thermodynamic database.
    shape_factor : float
        Volumetric shape factor of the solid crystals [-].

    Returns
    -------
    Cake
        Cake with real liquid and solid phase collaborators attached.
    """
    x_distrib = np.linspace(10.0, 500.0, 40)  # [um]
    mean_size = 200.0  # [um]
    size_std = 60.0  # [um]
    # Peak number density chosen in #158 to produce a litre-scale cake [#/um].
    peak_density = 8.0e5  # [#/um]
    distrib = peak_density * np.exp(
        -0.5 * ((x_distrib - mean_size) / size_std)**2)  # [#/um]

    liquid = LiquidPhase(thermo_path, mass_frac=LIQUID_MASS_FRAC)
    solid = SolidPhase(
        thermo_path,
        mass_frac=SOLID_MASS_FRAC,
        x_distrib=x_distrib,
        distrib=distrib,
        kv=shape_factor,
    )

    cake = Cake()
    cake.Phases = (liquid, solid)
    return cake


def test_cake_alpha_uses_phase_shape_factor_and_remains_invariant(thermo_path):
    """Verify the phase-owned scalar factor cancels from ``alpha``.

    Parameters
    ----------
    thermo_path : str
        Path to the pure-component thermodynamic database.
    """
    cake = _make_cake(thermo_path, SHAPE_FACTOR_PROBES[0])
    alpha_values = []  # [m/kg]
    micrometer_to_meter = 1e-6  # [m/um], exact unit conversion
    size_grid = cake.Solid_1.x_distrib * micrometer_to_meter  # [m]
    node_sizes = (size_grid[:-1] + size_grid[1:]) / 2  # [m]
    size_cubes = node_sizes**3  # [m**3]

    for shape_factor in SHAPE_FACTOR_PROBES:
        observed_shape_factor = MagicMock()
        observed_shape_factor.__mul__.return_value = (
            shape_factor * size_cubes
        )  # [m**3]
        with patch.object(
            type(cake.Solid_1),
            "kv",
            new_callable=PropertyMock,
            create=True,
            return_value=observed_shape_factor,
        ) as phase_shape_factor:
            alpha_values.append(cake.get_alpha())

        phase_shape_factor.assert_called_once_with()
        observed_shape_factor.__mul__.assert_called_once()
        multiplied_size_cubes = (
            observed_shape_factor.__mul__.call_args.args[0]
        )  # [m**3]
        np.testing.assert_array_equal(multiplied_size_cubes, size_cubes)

    alpha_values = np.asarray(alpha_values)  # [m/kg]
    relative_tolerance = np.finfo(float).eps  # [-], one binary64 epsilon
    assert alpha_values[0] == pytest.approx(
        EXPECTED_LEGACY_ALPHA,
        rel=relative_tolerance,
    )
    np.testing.assert_allclose(
        alpha_values,
        alpha_values[0],
        rtol=relative_tolerance,
        atol=0.0,
    )
