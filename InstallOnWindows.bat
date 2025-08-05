@echo OFF
echo -------------------------------------------------- 
echo Welcome to the PharmaPy installation wizard.
echo -------------------------------------------------- 
set /p env_name=Enter environment name: 
echo ------------------------------
echo Creating conda environment...
echo ------------------------------
call conda create -n %env_name% python=3.9 -c conda-forge -y
echo ------------------------------
echo Activating environment...
echo ------------------------------
call conda activate %env_name%
echo ------------------------------
echo Installing core dependencies...
echo ------------------------------
call conda install -c conda-forge numpy scipy matplotlib pandas cython -y
echo ------------------------------
echo Installing system dependencies for assimulo...
echo --------------------------------
echo Installing core dependencies...
echo --------------------------------
call pip install -r requirements.txt
echo ------------------------------
call conda install -c conda-forge sundials=5.8.0 superlu=5.2.2 openblas -y
echo ------------------------------
echo Attempting to install assimulo...
echo ------------------------------
call pip install -r requirements-assimulo.txt || echo "Assimulo installation failed - PharmaPy will work with limited functionality"
echo ----------------------
echo Installing PharmaPy...
echo ----------------------
REM Modern installation method (preferred)
call pip install -e .
REM Fallback to legacy method if pip install fails
if %ERRORLEVEL% neq 0 (
    echo Pip install failed, trying legacy method...
    call python setup.py develop
)
echo ------------------------------
echo Testing installation...
echo ------------------------------
python -c "import PharmaPy; print('✅ PharmaPy installed successfully')"
python -c "from PharmaPy import Utilities; print('✅ Core modules working')"
python -c "from PharmaPy import Reactors; print('✅ Simulation modules working')" || echo "⚠️ Simulation modules require assimulo"
echo ------------------------------
echo Installation complete!
echo ------------------------------
echo To activate this environment in the future, use:
echo conda activate %env_name%
echo.
echo To test the installation, run:
echo cd tests
echo python reactor_tests.py
echo.
pause
