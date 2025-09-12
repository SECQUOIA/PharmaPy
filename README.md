# PharmaPy

<img align="left" src="./doc/online_docs/images/PharmaPy_logo.jpeg" alt="PharmaPy_logo" height="250">

<!-- BEGIN Status badges -->
![GitHub all releases](https://img.shields.io/github/downloads/CryPTSys/PharmaPy/total)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.compchemeng.2021.107408-blue)](https://www.sciencedirect.com/science/article/abs/pii/S0098135421001861)
<!-- END Status badges -->

PharmaPy is a pythonic library for the analysis of pharmaceutical manufacturing systems.

It allows to simulate the dynamics of standalone, drug substance unit operations in a variety of operating modes (batch, continuous, semibatch). Also, PharmaPy facilitates setting up and simulating pharmaceutical **flowsheets**, i.e., interconnected unit operations running in one or more operation modes, offering flexibility to simulate end-to-end batch, end-to-end continuous, and hybrid operation schemes (combination of batch and/or continuous and semicontinuous unit operations).

<br clear="left"/>

## Getting started

### Quick Installation
For most users, the simplest installation method is:

```bash
# Create a conda environment (recommended)
conda create -n pharmapy python=3.9
conda activate pharmapy

# Install PharmaPy
pip install -e .
```

### Installation Options

#### Option 1: Using installation scripts (Automated)
Download and unzip the code from the release section, then:
- **Windows**: Run `InstallOnWindows.bat`
- **Linux/Mac**: Run `chmod +x InstallOnMac.sh && ./InstallOnMac.sh`

#### Option 2: Manual installation
Follow the detailed instructions in the `install_instructions.txt` file.

#### Option 3: Development installation
For developers:
```bash
# Clone the repository
git clone https://github.com/SECQUOIA/PharmaPy.git
cd PharmaPy

# Install with development dependencies
pip install -e ".[dev]"
```

### Testing the Installation
```bash
# Quick test
cd tests
python reactor_tests.py

# Comprehensive test suite
python run_tests.py

# Windows users can also use:
run_tests.bat

# Or with make:
make test          # Full test suite
make test-quick    # Basic functionality only
```

### Development and Testing
PharmaPy includes a comprehensive testing infrastructure:
- **Cross-platform testing**: Windows, macOS, Linux
- **Multi-Python support**: Python 3.9-3.12
- **CI/CD pipeline**: Automated testing via GitHub Actions
- **Multiple test runners**: Python scripts, batch files, Makefile targets

For developers:
```bash
# Install with development tools
pip install -e ".[dev]"

# Run code quality checks
make lint          # Code linting
make format        # Code formatting
```

Read our [documentation](https://pharmapy.readthedocs.io/en/latest/) or chat with the [PharmaPy Simulation Assistant](https://chatgpt.com/g/g-679bb3b5c5188191b26680b147a4f4a2-pharmapy-simulation-assistant) for more information on how to install and how to use PharmaPy.



