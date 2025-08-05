"""
Pytest configuration and shared fixtures for PharmaPy tests.
"""
import pytest
import os
import sys
import warnings
from pathlib import Path

# Add the PharmaPy package to the Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import PharmaPy modules
try:
    import PharmaPy
    from PharmaPy import Streams, Phases, Kinetics, Utilities
    PHARMAPY_AVAILABLE = True
except ImportError:
    PHARMAPY_AVAILABLE = False

# Check for optional dependencies
try:
    import assimulo
    ASSIMULO_AVAILABLE = True
    ASSIMULO_VERSION = assimulo.__version__
except ImportError:
    ASSIMULO_AVAILABLE = False
    ASSIMULO_VERSION = None

# Suppress warnings during tests
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


@pytest.fixture(scope="session", autouse=True)
def setup_test_environment():
    """Set up the test environment once per session."""
    # Ensure data directory exists
    data_dir = Path(__file__).parent.parent / "data"
    if not data_dir.exists():
        pytest.skip("Data directory not found")
    
    yield
    
    # Cleanup after all tests
    # Remove any temporary files created during testing
    pass


@pytest.fixture
def pharmapy_available():
    """Skip test if PharmaPy is not available."""
    if not PHARMAPY_AVAILABLE:
        pytest.skip("PharmaPy not available")
    return True


@pytest.fixture
def assimulo_available():
    """Skip test if assimulo is not available."""
    if not ASSIMULO_AVAILABLE:
        pytest.skip("assimulo not available - simulation tests require assimulo")
    return True


@pytest.fixture
def data_path():
    """Return the path to the test data directory."""
    return Path(__file__).parent.parent / "data"


@pytest.fixture
def temp_dir(tmp_path):
    """Create a temporary directory for test outputs."""
    return tmp_path


@pytest.fixture
def sample_liquid_stream():
    """Create a sample liquid stream for testing."""
    if not PHARMAPY_AVAILABLE:
        pytest.skip("PharmaPy not available")
    
    from PharmaPy.Streams import LiquidStream
    from PharmaPy.Phases import LiquidPhase
    
    # Create a simple test stream
    phase = LiquidPhase(["water", "ethanol"], [0.7, 0.3])
    stream = LiquidStream(phase, flow_rate=100.0, temperature=298.15)
    return stream


@pytest.fixture
def sample_reaction_kinetics():
    """Create sample reaction kinetics for testing."""
    if not PHARMAPY_AVAILABLE:
        pytest.skip("PharmaPy not available")
    
    from PharmaPy.Kinetics import RxnKinetics
    
    # Create simple kinetics (A -> B)
    kinetics = RxnKinetics(
        components=["A", "B"],
        reactions=["A -> B"],
        rate_constants=[1.0]
    )
    return kinetics


@pytest.fixture(scope="session")
def test_data_files():
    """List available test data files."""
    data_dir = Path(__file__).parent.parent / "data"
    if not data_dir.exists():
        return []
    
    json_files = list(data_dir.glob("*.json"))
    csv_files = list(data_dir.glob("*.csv"))
    
    return {
        "json": json_files,
        "csv": csv_files,
        "all": json_files + csv_files
    }


def pytest_configure(config):
    """Configure pytest with custom markers and settings."""
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as integration test"
    )
    config.addinivalue_line(
        "markers", "unit: mark test as unit test"
    )
    config.addinivalue_line(
        "markers", "assimulo: mark test as requiring assimulo"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers automatically."""
    for item in items:
        # Add slow marker to tests that might be slow
        if "integration" in item.keywords or "flowsheet" in item.keywords:
            item.add_marker(pytest.mark.slow)
        
        # Add assimulo marker to tests that use assimulo
        if "assimulo" in str(item.fspath).lower() or any(
            "assimulo" in str(arg) for arg in item.fixturenames
        ):
            item.add_marker(pytest.mark.assimulo)


def pytest_report_header(config):
    """Add custom header to pytest report."""
    header_lines = [
        f"PharmaPy Available: {PHARMAPY_AVAILABLE}",
        f"Assimulo Available: {ASSIMULO_AVAILABLE}",
    ]
    
    if ASSIMULO_AVAILABLE:
        header_lines.append(f"Assimulo Version: {ASSIMULO_VERSION}")
    
    return header_lines


def pytest_runtest_setup(item):
    """Setup for each test item."""
    # Skip assimulo tests if assimulo is not available
    if item.get_closest_marker("assimulo") and not ASSIMULO_AVAILABLE:
        pytest.skip("assimulo not available")
