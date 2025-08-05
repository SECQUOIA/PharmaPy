# PharmaPy Installation Guide

PharmaPy supports multiple installation methods. For the best compatibility, especially in CI/CD environments, we recommend using conda-forge.

## Quick Start (Recommended)

### Method 1: Complete Environment Setup
```bash
conda env create -f environment.yml
conda activate pharmapy
```

### Method 2: Conda-forge Installation
```bash
conda install -c conda-forge numpy scipy matplotlib pandas cython assimulo
```

### Method 3: Using Requirements Files
```bash
# Core dependencies only (pip)
pip install -r requirements.txt

# All dependencies including assimulo (conda)
conda install -c conda-forge --file requirements-conda.txt
```

## Installation Methods Comparison

| Method | Use Case | Assimulo Support | CI/CD Ready |
|--------|----------|------------------|-------------|
| `environment.yml` | New development setup | ✅ Full | ✅ Yes |
| `requirements-conda.txt` | Conda environments | ✅ Full | ✅ Yes |
| `requirements.txt` | pip-only environments | ❌ Manual | ⚠️ Limited |

## For GitHub Actions

Use the conda-forge approach for consistent builds:

```yaml
- uses: conda-incubator/setup-miniconda@v2
  with:
    miniforge-version: latest
    activate-environment: pharmapy
    environment-file: environment.yml
```

Or install packages directly:

```yaml
- uses: conda-incubator/setup-miniconda@v2
- run: conda install -c conda-forge --file requirements-conda.txt
```

## Assimulo (Simulation Engine)

Assimulo is required for PharmaPy's simulation features but has complex dependencies:

- ✅ **Available on conda-forge** (recommended)
- ❌ **NOT reliable on PyPI** (incomplete source distributions)
- 🔧 **Requires Fortran libraries** (pre-compiled in conda-forge)

### Working Versions
- assimulo: 3.7.1+ (conda-forge)
- numpy: 2.0+ (conda-forge)
- sundials: 7.4.0+ (auto-installed)

## Troubleshooting

### Common Issues
1. **Assimulo import errors**: Use conda-forge, not pip
2. **Dependency conflicts**: Use `environment.yml` for clean setup
3. **CI build failures**: Ensure all packages from conda-forge

### Verification
```bash
python validate_local.py
```

This runs the PharmaPy validation suite to confirm all features work correctly.
