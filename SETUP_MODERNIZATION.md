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
# Simple installation
pip install -e .

# Or using the automated scripts
# Windows: .\InstallOnWindows.bat
# Linux/Mac: ./InstallOnMac.sh
```

### For Developers
```bash
# With development tools
pip install -e ".[dev]"

# Or specific dependency groups
pip install -e ".[test]"    # Testing tools
pip install -e ".[docs]"    # Documentation tools
```

### Legacy Method (still supported)
```bash
python setup.py develop
```

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

## CI/CD Pipeline

The GitHub Actions workflow automatically:
- Tests on multiple Python versions (3.9-3.12)
- Tests on multiple operating systems (Ubuntu, Windows, macOS)
- Builds distribution packages
- Validates package integrity

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
