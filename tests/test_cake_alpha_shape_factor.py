"""Regression tests for the shape factor used by ``Cake.get_alpha``.

Issue #160 records that the scalar volumetric shape factor cancels from the
normalized crystal-volume weights. The hydraulic resistance ``alpha`` should
therefore be independent of any positive scalar ``kv``, while the method must
still obtain that factor from its attached solid phase rather than a literal.
"""

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

# Four binary64 epsilon allow the weighted sum to reassociate across NumPy
# builds while remaining a few units in the last place around the reference.
ALPHA_RELATIVE_TOLERANCE = 4 * np.finfo(float).eps  # [-]

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


def test_cake_alpha_is_invariant_to_real_phase_shape_factor(thermo_path):
    """Verify real solid phases give the pinned, shape-invariant ``alpha``.

    Parameters
    ----------
    thermo_path : str
        Path to the pure-component thermodynamic database.
    """
    alpha_values = np.array(
        [
            _make_cake(thermo_path, shape_factor).get_alpha()
            for shape_factor in SHAPE_FACTOR_PROBES
        ]
    )  # [m/kg]
    assert alpha_values[0] == pytest.approx(
        EXPECTED_LEGACY_ALPHA,
        rel=ALPHA_RELATIVE_TOLERANCE,
    )
    np.testing.assert_allclose(
        alpha_values,
        EXPECTED_LEGACY_ALPHA,
        rtol=ALPHA_RELATIVE_TOLERANCE,
        atol=0.0,
    )


@pytest.mark.parametrize(
    "invalid_shape_factor",
    (0.0, -0.5, np.nan, np.inf),
)
def test_cake_alpha_rejects_invalid_solid_phase_shape_factor(
    thermo_path,
    invalid_shape_factor,
):
    """Verify ``alpha`` validates ``kv`` from the attached solid phase.

    Parameters
    ----------
    thermo_path : str
        Path to the pure-component thermodynamic database.
    invalid_shape_factor : float
        Nonpositive or nonfinite volumetric shape factor [-].
    """
    cake = _make_cake(thermo_path, LEGACY_SHAPE_FACTOR)
    cake.Solid_1.kv = invalid_shape_factor  # [-]

    with pytest.raises(
        ValueError,
        match="shape factor 'kv' must be a finite, positive scalar",
    ):
        cake.get_alpha()
