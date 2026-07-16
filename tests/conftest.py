import importlib
from pathlib import Path
import sys
from types import ModuleType

import pytest


TESTS_ROOT = Path(__file__).resolve().parent
REPO_ROOT = TESTS_ROOT.parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _has_assimulo():
    try:
        import assimulo  # noqa: F401
    except ImportError:
        return False
    return True


def _stub_assimulo_modules(monkeypatch, *, solvers, problem):
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
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        if exc.name != "assimulo":
            raise

        _stub_assimulo_modules(monkeypatch, solvers=solvers, problem=problem)
        return importlib.import_module(module_name)


@pytest.fixture(scope="session")
def data_path():
    return {
        "integration": TESTS_ROOT / "integration" / "data",
        "flowsheet": TESTS_ROOT / "Flowsheet" / "data",
    }


def pytest_collection_modifyitems(config, items):
    if _has_assimulo():
        return

    skip_assimulo = pytest.mark.skip(
        reason="assimulo is not installed; solver-backed integration tests skipped"
    )
    for item in items:
        if "assimulo" in item.keywords:
            item.add_marker(skip_assimulo)
