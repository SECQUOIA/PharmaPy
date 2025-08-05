#!/bin/bash
echo -------------------------------------------------- 
echo Welcome to the PharmaPy installation wizard.
echo -------------------------------------------------- 
echo Enter environment name:
read env_name 
echo ------------------------------
echo Creating conda environment...
echo ------------------------------
conda create -n $env_name python=3.9 --file requirements.txt -c conda-forge -y
echo ------------------------------
echo Activating conda environment...
echo ------------------------------
conda init
eval "$(conda shell.bash hook)"
conda activate $env_name
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
echo Installation complete!
echo ------------------------------
echo "To activate this environment in the future, use:"
echo "conda activate $env_name"
echo
echo "To test the installation, run:"
echo "cd tests"
echo "python reactor_tests.py"
