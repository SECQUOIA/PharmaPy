# PharmaPy Setup Pipeline Modernization

## Summary of Changes

This document outlines the modernization of PharmaPy's setup pipeline while maintaining backward compatibility and Windows support.

## Files Modified/Created

### Core Build Configuration
- **`pyproject.toml`** (new): Modern Python packaging standard (PEP 517/518)
- **`setup.py`** (updated): Enhanced for backward compatibility while recommending pyproject.toml
- **`setup.cfg`** (updated): Declarative configuration file
- **`MANIFEST.in`** (updated): Ensures all necessary files are included in distributions

### Installation Scripts
- **`InstallOnWindows.bat`** (updated): 
  - Uses `pip install -e .` as primary method
  - Falls back to `python setup.py develop` if pip fails
  - Added confirmation messages and better error handling
  
- **`InstallOnMac.sh`** (updated):
  - Uses `pip install -e .` as primary method  
  - Falls back to `python setup.py develop` if pip fails
  - Improved environment activation and error handling

### Documentation
- **`install_instructions.txt`** (updated): Added modern installation methods
- **`README.md`** (updated): Clear installation options for different user types
- **`requirements-dev.txt`** (new): Development dependencies

### Development Tools  
- **`Makefile`** (updated): Common development tasks for Unix-like systems
- **`make.bat`** (updated): Windows equivalent of Makefile
- **`.github/workflows/test.yml`** (new): GitHub Actions CI/CD pipeline

## Key Features

### 1. Modern Standards Compliance
- Uses PEP 517/518 build system
- SPDX license specification  
- Modern dependency specification

### 2. Backward Compatibility
- `setup.py` still works for legacy workflows
- Installation scripts support both old and new methods
- Graceful fallback mechanisms

### 3. Cross-Platform Support
- Windows batch files with error handling
- Unix shell scripts with proper environment handling
- Platform-specific instructions in documentation

### 4. Development Experience
- Optional dependency groups (`[dev]`, `[test]`, `[docs]`)
- Code formatting and linting tools integration
- Automated testing pipeline

## Installation Methods

### For End Users
```bash
# Basic installation (core functionality)
pip install -e .

# With simulation capabilities (requires compatible assimulo)
pip install -e ".[simulation]"

# Or using the automated scripts
# Windows: .\InstallOnWindows.bat
# Linux/Mac: ./InstallOnMac.sh
```

### For Developers
```bash
# With development tools
pip install -e ".[dev]"

# With all dependencies (including simulation)
pip install -e ".[all]"

# Or specific dependency groups
pip install -e ".[test]"    # Testing tools
pip install -e ".[docs]"    # Documentation tools
```

### Legacy Method (still supported)
```bash
python setup.py develop
```

## Dependency Management

PharmaPy uses a tiered dependency approach:

### Core Dependencies (always installed)
- numpy: Numerical computing
- scipy: Scientific computing  
- matplotlib: Plotting
- pandas: Data manipulation
- cython: Python extensions

### Optional Dependencies
- **simulation**: `assimulo` for ODE solving (may have compatibility issues with newer Python versions)
- **dev**: Development tools (black, flake8, isort, etc.)
- **test**: Testing frameworks (pytest, pytest-cov)
- **docs**: Documentation tools (sphinx, themes)

### Compatibility Notes
- `assimulo` may not be compatible with Python 3.10+ or recent Cython versions
- For simulation features requiring assimulo, consider using Python 3.9 with older Cython versions
- Core PharmaPy functionality works without assimulo

## Build and Distribution

### Building the Package
```bash
# Install build tools
pip install build

# Build distribution packages
python -m build
```

### Using Development Tools
```bash
# Unix/Linux/Mac
make install-dev    # Install with dev dependencies
make test          # Run tests
make build         # Build package
make clean         # Clean artifacts

# Windows
make.bat install-dev
make.bat test
make.bat build
make.bat clean
```

## Testing Infrastructure

### GitHub Actions CI/CD Pipeline

#### Comprehensive Testing Workflows
1. **`test.yml`**: Main CI/CD pipeline with multiple test jobs
   - **Installation Tests**: Tests different installation methods across OS/Python versions
   - **Unit Tests**: Runs functional tests (limited to Python 3.9 for assimulo compatibility)
   - **Script Tests**: Tests automated installation scripts with Conda
   - **Build Tests**: Tests package building and distribution
   - **Code Quality**: Runs linting and formatting checks

2. **`installation-tests.yml`**: Dedicated installation script testing
   - Tests both Windows and Unix installation scripts
   - Tests manual installation procedures  
   - Validates functionality after installation

#### Test Matrix Coverage
- **Operating Systems**: Ubuntu, Windows, macOS
- **Python Versions**: 3.9-3.12 (with 3.9 focus for full functionality)
- **Installation Methods**: pip (basic/dev), setup.py develop
- **Dependency Scenarios**: Core dependencies vs. full simulation stack

### Local Testing Tools

#### Test Runners
1. **`run_tests.py`**: Comprehensive Python test runner
   - Graceful handling of optional dependencies
   - Modular test execution (imports, installation, reactor, flowsheet, pytest)
   - Cross-platform compatibility
   - Detailed reporting and timing

2. **`run_tests.bat`**: Windows batch test runner
   - Native Windows command prompt support
   - Same test coverage as Python runner
   - Visual status indicators

3. **`test_discovery.py`**: Advanced test discovery tool
   - Finds and runs unittest-based tests
   - Handles direct execution tests
   - Module import testing

#### Makefile Integration
- **Unix/Linux/Mac**: `make test`, `make test-quick`, `make test-reactor`, etc.
- **Windows**: `make.bat test`, `make.bat test-quick`, etc.
- Multiple test granularities (quick, full, specific modules)

### Test Configuration

#### Pytest Configuration (`pytest.ini`)
- Structured test discovery patterns
- Markers for different test types (unit, integration, slow)
- Logging and output configuration
- Cross-platform compatibility settings

#### Dependency Handling
- **Core Tests**: Work with basic dependencies (numpy, scipy, matplotlib, pandas)
- **Simulation Tests**: Require assimulo (Python 3.9 recommended)
- **Optional Features**: Graceful degradation when dependencies unavailable
- **Development Tools**: Additional linting and formatting tools

### Running Tests

#### Quick Testing
```bash
# Basic functionality test
make test-quick
# or
python run_tests.py --skip-flowsheet --skip-pytest

# Windows
make.bat test-quick
# or
run_tests.bat
```

#### Comprehensive Testing
```bash
# Full test suite
make test
# or  
python run_tests.py

# Specific test categories
make test-reactor      # Reactor functionality
make test-flowsheet    # Flowsheet functionality
make test-imports      # Import testing
```

#### CI/CD Integration
- Automated testing on every push/PR
- Multiple environment validation
- Installation method verification
- Dependency compatibility checking

### Test Coverage

#### Current Test Status
- ✓ **Package Installation**: Multiple methods validated
- ✓ **Core Imports**: Basic PharmaPy functionality
- ! **Simulation Features**: Require assimulo (compatibility issues)
- ✓ **Cross-Platform**: Windows, macOS, Linux support
- ✓ **Multi-Python**: 3.9-3.12 installation testing
- ✓ **Build System**: Modern packaging standards

#### Known Limitations
- `assimulo` compatibility issues with Python 3.10+
- Some tests require external solvers
- Simulation tests may need specific environment setup

## Benefits

1. **Future-Proof**: Uses modern Python packaging standards
2. **Backward Compatible**: Existing workflows continue to work
3. **Cross-Platform**: Robust Windows, Mac, and Linux support
4. **Developer Friendly**: Enhanced development tools and workflows
5. **Automated**: CI/CD pipeline for quality assurance
6. **Flexible**: Multiple installation methods for different use cases

## Migration Path

- **Immediate**: All current installation methods continue to work
- **Recommended**: Start using `pip install -e .` for new installations
- **Future**: Gradually migrate to pyproject.toml-based workflows
- **Legacy**: setup.py remains available for compatibility

This modernization ensures PharmaPy can leverage modern Python tooling while maintaining compatibility with existing workflows and providing excellent cross-platform support.
