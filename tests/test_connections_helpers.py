from types import SimpleNamespace

import numpy as np
import pytest

from PharmaPy.Connections import (
    convert_str_flowsheet,
    get_input_dict,
    get_missing_field,
    get_remaining_states,
    interpolate_inputs,
)


pytestmark = pytest.mark.unit


def test_convert_str_flowsheet_builds_adjacency_map_from_sequence():
    flowsheet = "Feed --> Reactor --> Filter"

    result = convert_str_flowsheet(flowsheet)

    assert result == {
        "Feed": ["Reactor"],
        "Reactor": ["Filter"],
        "Filter": [],
    }


@pytest.mark.parametrize(
    ("dim_state", "n_times", "expected"),
    [
        (1, 1, 0),
        (3, 1, np.zeros(3)),
        (1, 2, np.zeros(2)),
        (2, 3, np.zeros((3, 2))),
    ],
)
def test_get_missing_field_returns_zero_with_expected_shape(
    dim_state,
    n_times,
    expected,
):
    result = get_missing_field(dim_state, n_times)

    np.testing.assert_allclose(result, expected)


def test_get_input_dict_splits_vector_input_by_state_dimensions():
    input_data = np.array([1.0, 2.0, 300.0])
    name_dict = {"Inlet": {"mole_conc": 2, "temp": 1}}

    result = get_input_dict(input_data, name_dict)

    np.testing.assert_allclose(result["Inlet"]["mole_conc"], [1.0, 2.0])
    assert result["Inlet"]["temp"] == pytest.approx(300.0)


def test_get_input_dict_splits_time_series_columns_by_state_dimensions():
    input_data = np.array([
        [1.0, 2.0, 300.0],
        [3.0, 4.0, 310.0],
    ])
    name_dict = {"Inlet": {"mole_conc": 2, "temp": 1}}

    result = get_input_dict(input_data, name_dict)

    np.testing.assert_allclose(result["Inlet"]["mole_conc"], input_data[:, :2])
    np.testing.assert_allclose(result["Inlet"]["temp"], [300.0, 310.0])


def test_interpolate_inputs_evaluates_vector_times_without_solver_stack():
    inlet_time = np.array([0.0, 1.0, 2.0])
    inlet_values = np.array([
        [0.0, 0.0],
        [1.0, 10.0],
        [2.0, 20.0],
    ])

    result = interpolate_inputs(
        np.array([0.5, 1.5]),
        inlet_time,
        inlet_values,
    )

    np.testing.assert_allclose(result, [[0.5, 5.0], [1.5, 15.0]])


def test_get_remaining_states_uses_stream_values_and_zero_defaults():
    stream = SimpleNamespace(
        mole_conc=np.array([1.0, 2.0]),
        temp=300.0,
        Solid=SimpleNamespace(mu_n=None),
    )
    state_dimensions = {
        "Inlet": {"mole_conc": 2, "temp": 1, "mass": 1},
        "Solid": {"mu_n": 3},
    }
    existing_inlets = {"Inlet": {"temp": 305.0}, "Solid": {}}

    result = get_remaining_states(
        state_dimensions,
        stream,
        existing_inlets,
        time=np.array([0.0, 1.0]),
    )

    np.testing.assert_allclose(
        result["Inlet"]["mole_conc"],
        [[1.0, 2.0], [1.0, 2.0]],
    )
    np.testing.assert_allclose(result["Inlet"]["mass"], [0.0, 0.0])
    np.testing.assert_allclose(result["Solid"]["mu_n"], np.zeros((2, 3)))
