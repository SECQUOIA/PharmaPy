"""Repository policy checks for prohibited test substitutes."""

import ast
import hashlib
from pathlib import Path

import pytest


pytestmark = pytest.mark.unit

TESTS_ROOT = Path(__file__).resolve().parent
FORBIDDEN_MODULES = frozenset(
    {"_pytest.monkeypatch", "mock", "pytest_mock", "unittest.mock"}
)
FORBIDDEN_CONSTRUCTORS = frozenset(
    {
        "AsyncMock",
        "MagicMock",
        "Mock",
        "MonkeyPatch",
        "NonCallableMagicMock",
        "NonCallableMock",
        "PropertyMock",
    }
)
FORBIDDEN_FIXTURES = frozenset({"mocker", "monkeypatch"})

# Issue #202 removes these exact legacy files in focused PRs. The content
# digests form a ratchet: any edit to a file that still uses monkeypatch makes
# the repository check fail, while a migration that removes every violation
# passes without maintaining an exemption.
LEGACY_MONKEYPATCH_FILE_DIGESTS = {
    "tests/test_batch_cryst_concentration_jacobian.py": (
        "df111caf9baee2bd10a3f983ef8f90d902f2b678505328125fa30c62331db0f2"
    ),
    "tests/test_deliquoring_particle_size_units.py": (
        "e9dc7ce86a1361c3f4b486586678d7f0f801aff3ad7bceab01ee9eed913a12d6"
    ),
    "tests/test_drying_energy_rate_basis.py": (
        "92c65c66296a5a57ee6e06fc695bc6eccf215dc37bc597aab49726c1c4b2fb88"
    ),
    "tests/test_drying_gas_balance.py": (
        "1efc9c1940200e60b39ecbf8cd020f049c6ccc57d3f9af070964abda8cd22c38"
    ),
    "tests/test_drying_latent_heat_factor.py": (
        "0c50a7437f94113114901fc630f5b167d0ec00546bdda7027b64805d22d00fd7"
    ),
    "tests/test_drying_model.py": (
        "d4592e91a578bcfec961530a659d3490b84d6b376231775e96d91f2ad773c9c5"
    ),
    "tests/test_dynamic_collector_state_ordering.py": (
        "10af90d6a1a1c682d6f3ba6de5e034b4b55d4f0139d6b5b930ed4e173c37d58c"
    ),
    "tests/test_dynamic_extractor_default_k_fun.py": (
        "cf9a39a039806799756b6d3ea2b9353ad623afbd8b58cec02414feacbd8216d3"
    ),
    "tests/test_dynamic_extractor_stage_efficiency.py": (
        "8a666ab0e76c7423de36c1019f5449a1bd882240932f5790e755974575c2e71e"
    ),
    "tests/test_exception_handling.py": (
        "22617523a4774cd3dfc40400778915494b3c0fc932fe07f0c626868887287d9a"
    ),
    "tests/test_extractor_modes.py": (
        "657659bd3ed2c53eed18935e15e45586ce4e190a75d203c9e353b4ac7ad3e49b"
    ),
    "tests/test_metamodeler_codegen.py": (
        "958a52435715c302ea81b3d5c17ce4dc14fa9ae3ea9100f30774c4860a33f02d"
    ),
    "tests/test_optional_assimulo_imports.py": (
        "e65a16b56b0620f02a7e831f995c2e07b8f6b70b8caf8eb480257dc56d5f1480"
    ),
    "tests/test_paramestim_calibration_fit_predict.py": (
        "39da4b9c6ff5e81bf625cb6c7cf44b05f231d53ad364a1ef457c6e8826bca304"
    ),
    "tests/test_paramestim_gradient_assembly.py": (
        "b73418cd050a3abaa428caced65cb78c654f25801f970c808b30d92bfbce16c8"
    ),
    "tests/test_phases_vapor_latent_heat_shape.py": (
        "e0c96b1d4fea3726d5dde9ed95cd14d9b0f39fd42adc83da28eee18e0d9c72d2"
    ),
    "tests/test_reactors_energy.py": (
        "d636392346ebbef16de3bd6abf13f4f9aa6c6cbec1e5891466a4162a27bd0298"
    ),
    "tests/test_reversible_jacobian_reverse_term.py": (
        "e568cb5ac7d3ede1bf4429974303872a00b8276b1083f1de72513bc37fb619a8"
    ),
    "tests/test_simexec_routing.py": (
        "7abab90e232de7f27614c51ea4d424346254dd4a03682dbaacefcf9c3b1f5afa"
    ),
}


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


def _mock_policy_violations(
    test_file: Path, display_path: Path | None = None
) -> list[str]:
    """Find prohibited substitute APIs in one test file.

    Parameters
    ----------
    test_file : pathlib.Path
        Python test source to inspect.
    display_path : pathlib.Path, optional
        Path to show in diagnostics. Defaults to the repository-relative path.

    Returns
    -------
    list[str]
        Human-readable violations with repository-relative paths and lines.
    """
    source = test_file.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(test_file))
    relative_path = display_path or test_file.relative_to(TESTS_ROOT.parent)
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
            imports_pytest_monkeypatch = module_name == "pytest" and any(
                imported.name == "MonkeyPatch" for imported in node.names
            )
            if (
                _is_forbidden_module(module_name)
                or imports_unittest_mock
                or imports_pytest_monkeypatch
            ):
                violations.append(
                    f"{relative_path}:{node.lineno}: from {module_name} import"
                )
        elif isinstance(node, ast.arg) and node.arg in FORBIDDEN_FIXTURES:
            violations.append(
                f"{relative_path}:{node.lineno}: prohibited fixture '{node.arg}'"
            )
        elif (
            isinstance(node, ast.arg)
            and node.annotation is not None
            and ast.unparse(node.annotation).endswith("MonkeyPatch")
        ):
            violations.append(
                f"{relative_path}:{node.lineno}: pytest.MonkeyPatch annotation"
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


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("from unittest.mock import Mock\nvalue = Mock()\n", "unittest.mock"),
        (
            "def test_example(monkeypatch):\n"
            "    monkeypatch.setattr(object, 'name', 1)\n",
            "prohibited fixture 'monkeypatch'",
        ),
        ("from pytest import MonkeyPatch\n", "from pytest import"),
        (
            "def test_example(patcher: pytest.MonkeyPatch):\n    pass\n",
            "pytest.MonkeyPatch annotation",
        ),
    ],
)
def test_policy_detects_prohibited_substitute_apis(
    tmp_path: Path, source: str, expected: str
) -> None:
    """Reject representative mock and monkeypatch spellings structurally.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated pytest-provided directory for the source sample.
    source : str
        Python source containing a prohibited test-substitute API.
    expected : str
        Diagnostic fragment expected from the structural scan.
    """
    test_file = tmp_path / "test_example.py"
    test_file.write_text(source, encoding="utf-8")

    violations = _mock_policy_violations(
        test_file, display_path=Path("tests/test_example.py")
    )

    assert any(expected in violation for violation in violations)


def test_tests_do_not_use_prohibited_substitutes() -> None:
    """Keep new mock and monkeypatch APIs out of the PharmaPy test suite."""
    test_files = sorted(TESTS_ROOT.rglob("*.py"))
    assert test_files, "No Python tests were discovered for the policy check"

    violations = []
    for test_file in test_files:
        file_violations = _mock_policy_violations(test_file)
        if not file_violations:
            continue

        relative_path = test_file.relative_to(TESTS_ROOT.parent).as_posix()
        legacy_digest = LEGACY_MONKEYPATCH_FILE_DIGESTS.get(relative_path)
        normalized_source = test_file.read_text(encoding="utf-8")
        current_digest = hashlib.sha256(
            normalized_source.encode("utf-8")
        ).hexdigest()
        if current_digest == legacy_digest:
            continue

        violations.extend(file_violations)

    assert not violations, "Prohibited test substitutes found:\n" + "\n".join(
        violations
    )
