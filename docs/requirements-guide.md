# PharmaPy Requirements Structure

This document explains PharmaPy's requirements management approach using separate requirements files for different use cases.

## Requirements Files

### Core Requirements (`requirements.txt`)
Contains essential dependencies needed for basic PharmaPy functionality:
- `numpy`, `scipy`, `matplotlib`, `pandas` - Core numerical computing
- `cython` - Build dependencies

**Installation:**
```bash
pip install -r requirements.txt
```

### Development Requirements (`requirements-dev.txt`)
Contains tools for PharmaPy development and testing:
- `pytest`, `pytest-cov` - Testing framework
- `sphinx`, `sphinx-rtd-theme` - Documentation
- `black`, `flake8`, `isort` - Code formatting and linting
- `build`, `twine` - Package building and uploading

**Installation:**
```bash
pip install -r requirements-dev.txt
```

### Simulation Requirements (`requirements-assimulo.txt`)
Contains optional simulation dependencies for advanced modeling:
- `assimulo` - Dynamic simulation engine

**Prerequisites:**
```bash
# Install system dependencies first via conda
conda install -c conda-forge sundials=5.8.0 superlu=5.2.2 openblas
```

**Installation:**
```bash
pip install -r requirements-assimulo.txt
```

### Full Requirements (`requirements-full.txt`) 
Legacy file containing all dependencies including assimulo. Use the individual requirements files above for better dependency management.

## Installation Approaches

### 1. Core Installation (Recommended for most users)
```bash
pip install -r requirements.txt
pip install -e .
```

### 2. Development Installation
```bash
pip install -r requirements.txt
pip install -r requirements-dev.txt
pip install -e .
```

### 3. Full Installation with Simulation Support
```bash
# Install system dependencies
conda install -c conda-forge sundials=5.8.0 superlu=5.2.2 openblas

# Install Python dependencies
pip install -r requirements.txt
pip install -r requirements-assimulo.txt
pip install -e .
```

### 4. Using Make Commands
```bash
make install          # Core installation
make install-dev      # Development installation  
make install-assimulo # Full installation with simulation support
```

## CI/CD Integration

The GitHub Actions workflows use these requirements files:
- **Core functionality testing**: Uses `requirements.txt`
- **Development tools testing**: Uses `requirements-dev.txt`  
- **Simulation testing**: Uses `requirements.txt` + `requirements-assimulo.txt`

## Benefits of This Approach

1. **Modular Dependencies**: Users can install only what they need
2. **Faster CI/CD**: Core tests don't need heavy simulation dependencies
3. **Clear Separation**: Development tools separate from runtime dependencies
4. **Graceful Degradation**: PharmaPy works without optional dependencies
5. **Version Control**: Each requirement set can be versioned independently
6. **Documentation**: Clear installation paths for different use cases

## Compatibility Notes

- **assimulo**: Best compatibility with Python 3.9
- **System Dependencies**: Requires conda for simulation functionality
- **Cross-Platform**: All requirements files work on Windows, Linux, and macOS
