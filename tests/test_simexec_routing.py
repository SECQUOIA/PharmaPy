"""Regression tests for graph-defined flowsheet stream routing.

The cases solve real batch mixers and transfer real liquid phases through real
``Connection`` and ``SimulationResult`` collaborators.
"""

from pathlib import Path

import numpy as np
import pytest

from PharmaPy.Containers import Mixer
from PharmaPy.Phases import LiquidPhase
from PharmaPy.SimExec import SimulationExec


pytestmark = pytest.mark.unit

PATH_PHYS = str(
    Path(__file__).resolve().parent
    / "Flowsheet"
    / "data"
    / "compound_database.json"
)
FEED_A_MASS_FRACTION = np.array([0.70, 0.20, 0.05, 0.03, 0.02])  # [-]
FEED_B_MASS_FRACTION = np.array([0.10, 0.20, 0.30, 0.20, 0.20])  # [-]


def _add_real_feed(mixer, temperature):
    """Assign two real liquid phases to a source mixer.

    Parameters
    ----------
    mixer : Mixer
        Source unit that receives batch liquid phases.
    temperature : float
        Temperature of the first liquid phase [K]. The second phase is 1 K
        warmer to exercise the real energy balance.
    """
    mixer.Inlets = LiquidPhase(
        PATH_PHYS,
        mass=1.0,  # [kg]
        mass_frac=FEED_A_MASS_FRACTION,
        temp=temperature,
    )
    mixer.Inlets = LiquidPhase(
        PATH_PHYS,
        mass=2.0,  # [kg]
        mass_frac=FEED_B_MASS_FRACTION,
        temp=temperature + 1.0,  # [K]
    )


def _configured_flowsheet(graph):
    """Build a real mixer flowsheet and feed every graph source.

    Parameters
    ----------
    graph : dict of str to list of str
        Directed acyclic flowsheet adjacency mapping.

    Returns
    -------
    SimulationExec
        Executor containing real ``Mixer`` units and source phases.
    """
    flowsheet = SimulationExec(PATH_PHYS, flowsheet=graph)
    destinations = {
        destination
        for successors in graph.values()
        for destination in successors
    }
    for index, name in enumerate(graph):
        mixer = Mixer()
        setattr(flowsheet, name, mixer)
        if name not in destinations:
            source_temperature = 300.0 + 5.0 * index  # [K]
            _add_real_feed(mixer, source_temperature)

    return flowsheet


def _connection_edges(flowsheet):
    """Map real connection endpoints back to their graph unit names.

    Parameters
    ----------
    flowsheet : SimulationExec
        Solved flowsheet containing real connection objects.

    Returns
    -------
    list of tuple of str
        Source-destination names in connection creation order.
    """
    unit_names = {
        id(getattr(flowsheet, name)): name
        for name in flowsheet.graph
    }
    return [
        (
            unit_names[id(connection.source_uo)],
            unit_names[id(connection.destination_uo)],
        )
        for connection in flowsheet.connections.values()
    ]


def _assert_connections_follow_graph(graph, created):
    """Assert every real handoff is an edge in the supplied graph.

    Parameters
    ----------
    graph : dict of str to list of str
        Directed acyclic flowsheet adjacency mapping.
    created : sequence of tuple of str
        Observed source-destination handoffs.
    """
    misrouted = [
        (source, destination)
        for source, destination in created
        if destination not in graph[source]
    ]
    assert not misrouted


def test_stream_handoff_follows_graph_edges_not_execution_order():
    """Route real phase handoffs through graph edges, not adjacent units."""
    graph = {"A": ["C"], "B": ["C"], "C": ["D"], "D": []}
    flowsheet = _configured_flowsheet(graph)

    flowsheet.SolveFlowsheet(verbose=False)
    created = _connection_edges(flowsheet)

    _assert_connections_follow_graph(graph, created)
    assert set(created) == {("A", "C"), ("B", "C"), ("C", "D")}
    assert len(flowsheet.C.Inlets) == 2
    assert len(flowsheet.D.Inlets) == 1


def test_stream_handoff_supports_fanout_to_multiple_successors():
    """Copy a real source phase to every graph successor."""
    graph = {"A": ["B", "C"], "B": ["D"], "C": ["D"], "D": []}
    flowsheet = _configured_flowsheet(graph)

    flowsheet.SolveFlowsheet(verbose=False)
    created = _connection_edges(flowsheet)

    _assert_connections_follow_graph(graph, created)
    assert set(created) == {
        ("A", "B"),
        ("A", "C"),
        ("B", "D"),
        ("C", "D"),
    }
    assert len(flowsheet.B.Inlets) == 1
    assert len(flowsheet.C.Inlets) == 1
    assert len(flowsheet.D.Inlets) == 2


def test_already_solved_branch_honors_pick_units():
    """Filter a real pre-solved unit's handoffs to selected destinations."""
    graph = {"A": ["B", "C"], "B": ["D"], "C": ["D"], "D": []}
    flowsheet = _configured_flowsheet(graph)
    _add_real_feed(flowsheet.B, 307.0)  # [K]
    flowsheet.B.solve_unit()

    flowsheet.SolveFlowsheet(pick_units=["A", "C"], verbose=False)
    created = _connection_edges(flowsheet)

    _assert_connections_follow_graph(graph, created)
    assert created == [("A", "C")]
    assert len(flowsheet.C.Inlets) == 1
    assert flowsheet.D.Inlets == []
