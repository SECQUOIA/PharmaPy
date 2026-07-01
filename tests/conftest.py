from pathlib import Path
import sys

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
