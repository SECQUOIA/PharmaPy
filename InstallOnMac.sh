#!/bin/bash
echo -------------------------------------------------- 
echo Welcome to the PharmaPy installation wizard.
echo -------------------------------------------------- 
echo Enter environment name:
read env_name 
echo ------------------------------
echo Creating conda environment...
echo ------------------------------
conda create -n $env_name python=3.9 -c conda-forge -y
echo ------------------------------
echo Activating conda environment...
echo ------------------------------
conda init
eval "$(conda shell.bash hook)"
conda activate $env_name
echo ------------------------------
echo Installing core dependencies...
echo ------------------------------
conda install -c conda-forge numpy scipy matplotlib pandas cython -y
echo ------------------------------
echo Installing system dependencies for assimulo...
echo ------------------------------
# Install dependencies needed for assimulo
echo --------------------------------
echo Installing core dependencies...
echo --------------------------------
pip install -r requirements.txt

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux
    conda install -c conda-forge sundials=5.8.0 superlu=5.2.2 openblas -y
elif [[ "$OSTYPE" == "darwin"* ]]; then
    # Mac OSX
    conda install -c conda-forge sundials=5.8.0 superlu=5.2.2 openblas -y
fi
echo ------------------------------
echo Attempting to install assimulo...
echo ------------------------------
pip install -r requirements-assimulo.txt || echo "! Assimulo installation failed - PharmaPy will work with limited functionality"
echo ----------------------
echo Installing PharmaPy...
echo ----------------------
# Modern installation method (preferred)
pip install -e .
# Check if pip install was successful
if [ $? -ne 0 ]; then
    echo "Pip install failed, trying legacy method..."
    python setup.py develop
fi
echo ------------------------------
echo Testing installation...
echo ------------------------------
python -c "import PharmaPy; print('✓ PharmaPy installed successfully')"
python -c "from PharmaPy import Utilities; print('✓ Core modules working')"
python -c "from PharmaPy import Reactors; print('✓ Simulation modules working')" || echo "! Simulation modules require assimulo"
echo ------------------------------
echo Installation complete!
echo ------------------------------
echo "To activate this environment in the future, use:"
echo "conda activate $env_name"
echo
echo "To test the installation, run:"
echo "cd tests"
echo "python reactor_tests.py"
