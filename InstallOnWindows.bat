@echo OFF
echo -------------------------------------------------- 
echo Welcome to the PharmaPy installation wizard.
echo -------------------------------------------------- 
set /p env_name=Enter environment name: 
echo ------------------------------
echo Creating conda environment...
echo ------------------------------
call conda create -n %env_name% python=3.9 --file requirements.txt -c conda-forge -y
echo ------------------------------
echo Activating environment...
echo ------------------------------
call conda activate %env_name%
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
