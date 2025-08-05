@echo off
setlocal

if "%1"=="" goto help
if "%1"=="help" goto help
if "%1"=="install" goto install
if "%1"=="install-dev" goto install-dev
if "%1"=="test" goto test
if "%1"=="test-quick" goto test-quick
if "%1"=="test-reactor" goto test-reactor
if "%1"=="test-flowsheet" goto test-flowsheet
if "%1"=="test-imports" goto test-imports
if "%1"=="test-all" goto test-all
if "%1"=="clean" goto clean
if "%1"=="build" goto build
if "%1"=="upload" goto upload
if "%1"=="docs" goto docs
if "%1"=="format" goto format
if "%1"=="lint" goto lint
goto invalid

:help
echo Available commands:
echo   install        Install the package in development mode
echo   install-dev    Install with development dependencies
echo   test           Run all tests
echo   test-quick     Run quick tests (imports + reactor tests)
echo   test-reactor   Run reactor tests only
echo   test-flowsheet Run flowsheet tests only
echo   test-imports   Test package imports only
echo   test-all       Run comprehensive test suite
echo   clean          Clean build artifacts
echo   build          Build the package
echo   upload         Upload to PyPI (requires credentials)
echo   docs           Build documentation
echo   format         Format code with black and isort
echo   lint           Run code linting
echo.
echo Usage: make.bat [command]
goto end

:install
pip install -e .
goto end

:install-dev
pip install -e ".[dev]"
goto end

:test
goto test-all

:test-quick
echo Running quick test suite...
python -c "import PharmaPy; print('✅ PharmaPy imports OK')"
if %ERRORLEVEL% neq 0 goto end
cd tests\integration
python reactor_tests.py
cd ..\..
goto end

:test-reactor
echo Running reactor tests...
cd tests\integration
python reactor_tests.py
cd ..\..
goto end

:test-flowsheet
echo Running flowsheet tests...
cd tests\Flowsheet
python flowsheet_tests.py
cd ..\..
goto end

:test-imports
echo Testing package imports...
python -c "import PharmaPy; print('✅ PharmaPy')" || echo "❌ PharmaPy import failed"
python -c "from PharmaPy import Utilities; print('✅ Utilities')" || echo "⚠️ Utilities import failed"
python -c "from PharmaPy import Reactors; print('✅ Reactors')" || echo "⚠️ Reactors import failed (may need assimulo)"
python -c "from PharmaPy import Streams; print('✅ Streams')" || echo "⚠️ Streams import failed (may need assimulo)"
python -c "from PharmaPy import Phases; print('✅ Phases')" || echo "⚠️ Phases import failed (may need assimulo)"
python -c "from PharmaPy import Kinetics; print('✅ Kinetics')" || echo "⚠️ Kinetics import failed (may need assimulo)"
goto end

:test-all
echo Running comprehensive test suite...
python run_tests.py
goto end

:clean
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
for /d %%i in (*.egg-info) do if exist "%%i" rmdir /s /q "%%i"
for /d /r . %%d in (__pycache__) do @if exist "%%d" rmdir /s /q "%%d"
del /s /q *.pyc 2>nul
if exist .pytest_cache rmdir /s /q .pytest_cache
if exist htmlcov rmdir /s /q htmlcov
if exist .coverage del .coverage
echo Cleaned build artifacts
goto end

:build
call :clean
python -m build
goto end

:upload
call :build
twine upload dist/*
goto end

:docs
cd doc
make html
cd ..
goto end

:format
black PharmaPy/
isort PharmaPy/
goto end

:lint
flake8 PharmaPy/
goto end

:invalid
echo Invalid command: %1
echo Use "make.bat help" to see available commands
goto end

:end
