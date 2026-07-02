import numpy as np
import pytest

from PharmaPy.Connections import topological_bfs
from PharmaPy.Interpolation import local_newton_interpolation
from PharmaPy.NameAnalysis import get_dict_states


def test_topological_bfs_orders_linear_graph():
    graph = {
        "feed": ["reactor"],
        "reactor": ["filter"],
        "filter": [],
    }

    in_degree, path = topological_bfs(graph)

    assert path == ["feed", "reactor", "filter"]
    assert all(count == 0 for count in in_degree.values())


def test_topological_bfs_leaves_cycle_with_incoming_edges():
    graph = {
        "a": ["b"],
        "b": ["a"],
    }

    in_degree, path = topological_bfs(graph)

    assert path == []
    assert in_degree == {"a": 1, "b": 1}


def test_get_dict_states_splits_composition_distribution_and_scalar():
    states = np.array([
        [1.0, 2.0, 10.0, 20.0, 30.0, 300.0],
        [3.0, 4.0, 40.0, 50.0, 60.0, 310.0],
    ])

    result = get_dict_states(
        ["mole_conc", "mu_n", "temp"],
        num_species=2,
        num_distr=3,
        states=states,
    )

    np.testing.assert_allclose(result["mole_conc"], states[:, :2])
    np.testing.assert_allclose(result["mu_n"], states[:, 2:5])
    np.testing.assert_allclose(result["temp"], states[:, 5])


def test_local_newton_interpolation_matches_quadratic():
    time = np.array([0.0, 1.0, 2.0, 3.0])
    values = time**2

    interpolated = local_newton_interpolation(1.5, time, values)

    assert interpolated == pytest.approx(2.25)
