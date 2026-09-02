"""Regression tests for DynamicExtractor default equilibrium routing."""

import numpy as np
import pytest

import PharmaPy.DynamicExtraction as dynamic_extraction
from PharmaPy.Phases import LiquidPhase
from PharmaPy.Streams import LiquidStream


pytestmark = pytest.mark.unit


def test_default_k_fun_uses_selected_activity_model(data_path):
    """The default DynamicExtractor k_fun returns gamma_light/gamma_heavy."""
    database_path = data_path["flowsheet"] / "compound_database.json"
    extractor = dynamic_extraction.DynamicExtractor(
        num_stages=2,
        gamma_model="ideal",
    )
    extractor.Phases = LiquidPhase(
        str(database_path),
        moles=10.0,  # [mol]
        mole_frac=np.array([0.20, 0.20, 0.20, 0.20, 0.20]),  # [-]
        temp=298.15,  # [K]
        pres=101325.0,  # [Pa]
    )

    x_light = np.array([
        [0.20, 0.10, 0.25, 0.15, 0.30],
        [0.10, 0.20, 0.15, 0.25, 0.30],
    ])  # [-]
    x_heavy = np.array([
        [0.15, 0.20, 0.10, 0.25, 0.30],
        [0.25, 0.10, 0.20, 0.15, 0.30],
    ])  # [-]
    temp = np.array([298.15, 301.15])  # [K]

    k_i = extractor.k_fun(x_light, x_heavy, temp)  # K_i [-]

    np.testing.assert_allclose(k_i, np.ones_like(x_light))
    assert extractor.gamma_model == "ideal"


def test_default_k_fun_rejects_unknown_activity_model():
    """DynamicExtractor rejects unknown default activity-model selectors."""
    # Lowercase ``uniquac`` is deliberate: activity-model selectors are
    # case-sensitive and must match the exact Phases.getActivityCoeff branch.
    with pytest.raises(ValueError, match="gamma_model must be one of"):
        dynamic_extraction.DynamicExtractor(
            num_stages=1,
            gamma_model="uniquac",
        )


def test_initialize_model_uses_default_equilibrium_with_real_phases(data_path):
    """Initialization solves the default equilibrium with real liquid objects."""
    database_path = data_path["flowsheet"] / "compound_database.json"
    # Representative extraction pair: the feed is A-rich and the solvent
    # stream is enriched in the database's designated solvent component.
    feed_mole_fraction = np.array([0.35, 0.20, 0.15, 0.15, 0.15])  # [-]
    solvent_mole_fraction = np.array([0.05, 0.20, 0.25, 0.25, 0.25])  # [-]
    holdup_mole_fraction = (
        feed_mole_fraction + solvent_mole_fraction
    ) / 2  # [-]
    extractor = dynamic_extraction.DynamicExtractor(
        num_stages=1,
        gamma_model="ideal",
    )
    extractor.Phases = LiquidPhase(
        str(database_path),
        moles=10.0,  # [mol]
        mole_frac=holdup_mole_fraction,
        temp=298.15,  # [K]
        pres=101325.0,  # [Pa]
    )
    extractor.Inlet = {
        "feed": LiquidStream(
            str(database_path),
            mole_flow=5.0,  # [mol/s]
            mole_frac=feed_mole_fraction,
            temp=298.15,  # [K]
            pres=101325.0,  # [Pa]
        ),
        "solvent": LiquidStream(
            str(database_path),
            mole_flow=3.0,  # [mol/s]
            mole_frac=solvent_mole_fraction,
            temp=298.15,  # [K]
            pres=101325.0,  # [Pa]
        ),
    }

    initial_states = extractor.initialize_model()

    np.testing.assert_allclose(initial_states["x_i"].sum(axis=1), 1.0)
    np.testing.assert_allclose(initial_states["y_i"].sum(axis=1), 1.0)
    np.testing.assert_allclose(initial_states["x_i"], initial_states["y_i"])
    assert extractor.fixed_vals["H_R"] > 0.0  # [mol]
    assert extractor.fixed_vals["H_E"] > 0.0  # [mol]
