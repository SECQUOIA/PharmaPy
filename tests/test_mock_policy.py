"""Repository policy checks for prohibited general-purpose mock objects."""

import ast
from pathlib import Path

import pytest


pytestmark = pytest.mark.unit

TESTS_ROOT = Path(__file__).resolve().parent
FORBIDDEN_MODULES = frozenset({"mock", "pytest_mock", "unittest.mock"})
FORBIDDEN_CONSTRUCTORS = frozenset(
    {
        "AsyncMock",
        "MagicMock",
        "Mock",
        "NonCallableMagicMock",
        "NonCallableMock",
        "PropertyMock",
    }
)


def _is_forbidden_module(module_name: str) -> bool:
    """Return whether an import belongs to a prohibited mock framework.

    Parameters
    ----------
    module_name : str
        Fully qualified imported module name.

    Returns
    -------
    bool
        ``True`` when the import is a prohibited module or submodule.
    """
    return any(
        module_name == forbidden
        or module_name.startswith(f"{forbidden}.")
        for forbidden in FORBIDDEN_MODULES
    )


def _mock_policy_violations(test_file: Path) -> list[str]:
    """Find prohibited mock imports and constructor calls in one test file.

    Parameters
    ----------
    test_file : pathlib.Path
        Python test source to inspect.

    Returns
    -------
    list[str]
        Human-readable violations with repository-relative paths and lines.
    """
    source = test_file.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(test_file))
    relative_path = test_file.relative_to(TESTS_ROOT.parent)
    violations = []

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for imported in node.names:
                if _is_forbidden_module(imported.name):
                    violations.append(
                        f"{relative_path}:{node.lineno}: import {imported.name}"
                    )
        elif isinstance(node, ast.ImportFrom):
            module_name = node.module or ""
            imports_unittest_mock = (
                module_name == "unittest"
                and any(imported.name == "mock" for imported in node.names)
            )
            if _is_forbidden_module(module_name) or imports_unittest_mock:
                violations.append(
                    f"{relative_path}:{node.lineno}: from {module_name} import"
                )
        elif isinstance(node, ast.arg) and node.arg == "mocker":
            violations.append(
                f"{relative_path}:{node.lineno}: pytest-mock fixture 'mocker'"
            )
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                constructor_name = node.func.id
            elif isinstance(node.func, ast.Attribute):
                constructor_name = node.func.attr
            else:
                continue

            if constructor_name in FORBIDDEN_CONSTRUCTORS:
                violations.append(
                    f"{relative_path}:{node.lineno}: {constructor_name}(...)"
                )

    return violations


def test_tests_do_not_use_general_purpose_mock_objects():
    """Keep prohibited mock frameworks out of the PharmaPy test suite."""
    test_files = sorted(TESTS_ROOT.rglob("*.py"))
    assert test_files, "No Python tests were discovered for the policy check"

    violations = [
        violation
        for test_file in test_files
        for violation in _mock_policy_violations(test_file)
    ]
    assert not violations, "Prohibited mock APIs found:\n" + "\n".join(violations)
