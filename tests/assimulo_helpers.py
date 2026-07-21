"""Helpers for tests that import optional Assimulo-backed modules."""

import importlib
import sys
from types import ModuleType


def _stub_assimulo_modules(monkeypatch, *, solvers, problem):
    """Register a minimal temporary Assimulo module tree.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest cleanup fixture used to restore ``sys.modules``.
    solvers : dict
        Solver attributes to expose on ``assimulo.solvers``.
    problem : dict
        Problem attributes to expose on ``assimulo.problem``.

    Notes
    -----
    This monkeypatch is limited to optional import availability in test
    environments where Assimulo is not installed. It does not implement solver
    behavior; solver-backed integration coverage remains in Assimulo-enabled
    environments.
    """
    assimulo = ModuleType("assimulo")

    solvers_module = ModuleType("assimulo.solvers")
    for name, value in solvers.items():
        setattr(solvers_module, name, value)

    problem_module = ModuleType("assimulo.problem")
    for name, value in problem.items():
        setattr(problem_module, name, value)

    exception = ModuleType("assimulo.exception")
    exception.TerminateSimulation = Exception

    monkeypatch.setitem(sys.modules, "assimulo", assimulo)
    monkeypatch.setitem(sys.modules, "assimulo.solvers", solvers_module)
    monkeypatch.setitem(sys.modules, "assimulo.problem", problem_module)
    monkeypatch.setitem(sys.modules, "assimulo.exception", exception)


def import_module_with_assimulo_stub(monkeypatch, module_name, *,
                                     solvers, problem):
    """Import a module, stubbing only optional Assimulo imports if needed.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest cleanup fixture used only when the import needs stubbing.
    module_name : str
        Dotted module path to import.
    solvers : dict
        Solver attributes to expose on ``assimulo.solvers``.
    problem : dict
        Problem attributes to expose on ``assimulo.problem``.

    Returns
    -------
    module
        Imported module.

    Raises
    ------
    ModuleNotFoundError
        If a dependency other than optional Assimulo is missing.

    Notes
    -----
    The temporary monkeypatch is justified because these unit tests exercise
    constructor and handoff logic that can run without Assimulo, while the full
    solver suite is only available when the optional Assimulo dependency is
    installed.
    """
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        if exc.name != "assimulo":
            raise

        _stub_assimulo_modules(monkeypatch, solvers=solvers, problem=problem)
        monkeypatch.delitem(sys.modules, module_name, raising=False)
        return importlib.import_module(module_name)
