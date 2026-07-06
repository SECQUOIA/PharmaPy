import os

import numpy as np
import pytest

import PharmaPy.SimExec as se


pytestmark = pytest.mark.unit

PATH_PHYS = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(se.__file__))),
    "tests",
    "Flowsheet",
    "data",
    "compound_database.json",
)


class _Result:
    time = np.array([0.0, 1.0])


class _StubUO:
    __module__ = "PharmaPy.Containers"

    def __init__(self, name):
        self._name = name
        self.outputs = {}
        self.result = _Result()

    def solve_unit(self, **kwargs):
        pass

    def flatten_states(self):
        pass


def _solve_and_record_connections(monkeypatch, graph):
    created = []

    class ConnRecorder:
        def __init__(self, source_uo, destination_uo):
            created.append((source_uo._name, destination_uo._name))

        def transfer_data(self):
            pass

    monkeypatch.setattr(se, "check_modeling_objects", lambda *args, **kwargs: None)
    monkeypatch.setattr(se, "SimulationResult", lambda *args, **kwargs: None)
    monkeypatch.setattr(se, "Connection", ConnRecorder)

    flst = se.SimulationExec(PATH_PHYS, flowsheet=graph)
    for name in graph:
        setattr(flst, name, _StubUO(name))

    flst.SolveFlowsheet(verbose=False)

    return created


def _assert_connections_follow_graph(graph, created):
    misrouted = [(source, dest) for source, dest in created
                 if dest not in graph[source]]

    assert not misrouted


def test_stream_handoff_follows_graph_edges_not_execution_order(monkeypatch):
    graph = {"A": ["C"], "B": ["C"], "C": ["D"], "D": []}

    created = _solve_and_record_connections(monkeypatch, graph)

    _assert_connections_follow_graph(graph, created)
    assert set(created) == {("A", "C"), ("B", "C"), ("C", "D")}


def test_stream_handoff_supports_fanout_to_multiple_successors(monkeypatch):
    graph = {"A": ["B", "C"], "B": ["D"], "C": ["D"], "D": []}

    created = _solve_and_record_connections(monkeypatch, graph)

    _assert_connections_follow_graph(graph, created)
    assert set(created) == {("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")}
