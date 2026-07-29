"""Import-boundary regressions for the optional Assimulo solver stack."""

import os
from pathlib import Path
import subprocess
import sys
import textwrap
from types import ModuleType

import pytest

import PharmaPy._assimulo as assimulo_backend

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[1]
AFFECTED_MODULES = (
    "PharmaPy.Containers",
    "PharmaPy.Crystallizers",
    "PharmaPy.Distillation",
    "PharmaPy.Drying_Model",
    "PharmaPy.DynamicExtraction",
    "PharmaPy.Evaporators",
    "PharmaPy.Reactors",
    "PharmaPy.SolidLiquidSep",
    "PharmaPy.ThreePhaseSettler",
)
LAZY_CONSTRUCTORS = (
    "CVode",
    "IDA",
    "Explicit_Problem",
    "Implicit_Problem",
)

IMPORT_BLOCKER = textwrap.dedent("""
    import builtins

    original_import = builtins.__import__

    def import_without_assimulo(
        name, globals=None, locals=None, fromlist=(), level=0
    ):
        if name == "assimulo" or name.startswith("assimulo."):
            raise ModuleNotFoundError(
                "Assimulo import blocked by regression test",
                name="assimulo",
            )
        return original_import(name, globals, locals, fromlist, level)

    builtins.__import__ = import_without_assimulo
    """)


def _run_without_assimulo(script, tmp_path):
    """Run Python source while rejecting imports of Assimulo.

    Parameters
    ----------
    script : str
        Python source to execute after installing the import blocker.
    tmp_path : pathlib.Path
        Temporary directory used for the Matplotlib configuration cache.

    Returns
    -------
    subprocess.CompletedProcess
        Completed child-process result with captured text output.
    """
    environment = os.environ.copy()
    environment["MPLCONFIGDIR"] = str(tmp_path)
    return subprocess.run(
        [sys.executable, "-c", f"{IMPORT_BLOCKER}\n{script}"],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def test_model_modules_import_without_assimulo(tmp_path):
    """Every solver-backed model module imports without the optional backend."""
    script = textwrap.dedent(f"""
        import importlib

        for module_name in {AFFECTED_MODULES!r}:
            importlib.import_module(module_name)
        """)
    result = _run_without_assimulo(script, tmp_path)

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("symbol_name", LAZY_CONSTRUCTORS)
def test_solver_construction_reports_missing_assimulo(symbol_name, tmp_path):
    """Lazy solver construction gives an actionable dependency error."""
    script = textwrap.dedent(f"""
        from PharmaPy._assimulo import {symbol_name}

        try:
            {symbol_name}(object())
        except ImportError as exc:
            expected = (
                "Assimulo is required for solver-backed PharmaPy simulations"
            )
            if expected not in str(exc):
                raise AssertionError(str(exc)) from exc
        else:
            raise AssertionError(
                "{symbol_name} construction unexpectedly succeeded"
            )
        """)
    result = _run_without_assimulo(script, tmp_path)

    assert result.returncode == 0, result.stderr


def test_missing_assimulo_symbol_reports_qualified_name(monkeypatch):
    """An incompatible install identifies the missing constructor."""
    stub_module = ModuleType("assimulo.solvers")

    def load_stub(module_name):
        assert module_name == "assimulo.solvers"
        return stub_module

    monkeypatch.setattr(assimulo_backend, "import_module", load_stub)

    with pytest.raises(
        ImportError,
        match=r"assimulo\.solvers\.CVode",
    ) as exc_info:
        assimulo_backend.CVode()

    assert isinstance(exc_info.value.__cause__, AttributeError)


def test_broken_assimulo_install_has_distinct_error(monkeypatch):
    """A backend loader failure is not misreported as a missing install."""
    loader_error = "SUNDIALS shared library could not be loaded"

    def fail_import(module_name):
        assert module_name == "assimulo.solvers"
        raise ImportError(loader_error)

    monkeypatch.setattr(assimulo_backend, "import_module", fail_import)

    with pytest.raises(
        ImportError,
        match="Assimulo is installed but could not be imported",
    ) as exc_info:
        assimulo_backend.CVode()

    assert loader_error in str(exc_info.value.__cause__)
