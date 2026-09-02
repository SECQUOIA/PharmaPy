"""Import-boundary regressions for the optional Assimulo solver stack."""

import os
from pathlib import Path
import subprocess
import sys
import textwrap

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

# Block Assimulo through a sys.meta_path finder rather than by patching
# builtins.__import__. importlib.import_module, which PharmaPy._assimulo uses to
# load the backend, resolves through importlib._bootstrap._gcd_import and never
# consults builtins.__import__, so a builtins patch lets the real Assimulo load
# whenever it is installed. A meta-path finder participates in that resolution
# and so blocks both statement imports and importlib.import_module.
#
# ModuleNotFoundError carries name="assimulo" because _load_assimulo_symbol
# distinguishes a missing install from a broken one on exactly that attribute.
IMPORT_BLOCKER = textwrap.dedent("""
    import sys

    class _AssimuloBlocker:
        def find_spec(self, name, path=None, target=None):
            if name == "assimulo" or name.startswith("assimulo."):
                raise ModuleNotFoundError(
                    "Assimulo import blocked by regression test",
                    name="assimulo",
                )
            return None

    sys.meta_path.insert(0, _AssimuloBlocker())
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
            # Both the missing-install and the incompatible-install branch embed
            # the message above, so assert the discriminator too: a missing
            # install is reported from ModuleNotFoundError, while an installed
            # backend lacking the symbol is reported from AttributeError and
            # adds "does not provide". Without this the test would still pass
            # while exercising the branch it is not named for.
            if not isinstance(exc.__cause__, ModuleNotFoundError):
                raise AssertionError(
                    "expected a missing-install error, got cause "
                    f"{{type(exc.__cause__).__name__}}: {{exc}}"
                ) from exc
        else:
            raise AssertionError(
                "{symbol_name} construction unexpectedly succeeded"
            )
        """)
    result = _run_without_assimulo(script, tmp_path)

    assert result.returncode == 0, result.stderr


@pytest.mark.assimulo
def test_missing_assimulo_symbol_reports_qualified_name():
    """The real backend identifies a requested constructor it lacks."""
    pytest.importorskip("assimulo.solvers")
    missing_symbol = "ConstructorThatAssimuloMustNotProvide"

    with pytest.raises(
        ImportError,
        match=rf"assimulo\.solvers\.{missing_symbol}",
    ) as exc_info:
        assimulo_backend._load_assimulo_symbol(
            "assimulo.solvers", missing_symbol
        )

    assert isinstance(exc_info.value.__cause__, AttributeError)


def test_broken_assimulo_install_has_distinct_error():
    """Loader failures classify as broken rather than missing installs."""
    loader_error = "SUNDIALS shared library could not be loaded"
    original_error = ImportError(loader_error)

    with pytest.raises(
        ImportError,
        match="Assimulo is installed but could not be imported",
    ) as exc_info:
        assimulo_backend._raise_assimulo_import_error(
            "assimulo.solvers", original_error
        )

    assert "assimulo.solvers" in str(exc_info.value)
    assert exc_info.value.__cause__ is original_error
