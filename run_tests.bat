@echo off
REM PharmaPy Test Runner for Windows
REM This script runs all PharmaPy tests

echo ================================================
echo PharmaPy Test Suite - Windows
echo ================================================

set FAILED_TESTS=0

echo.
echo Testing PharmaPy imports...
python -c "import PharmaPy; print('[OK] PharmaPy imported successfully')"
if %ERRORLEVEL% neq 0 (
    echo [FAIL] PharmaPy import failed
    set /a FAILED_TESTS+=1
)

echo.
echo Testing core module imports...
python -c "from PharmaPy import Reactors, Streams, Phases; print('[OK] Core modules imported successfully')"
if %ERRORLEVEL% neq 0 (
    echo [FAIL] Core module imports failed
    set /a FAILED_TESTS+=1
)

echo.
echo Running reactor tests...
cd tests\integration
python reactor_tests.py
if %ERRORLEVEL% neq 0 (
    echo [FAIL] Reactor tests failed
    set /a FAILED_TESTS+=1
) else (
    echo [OK] Reactor tests passed
)
cd ..\..

echo.
echo Running flowsheet tests...
cd tests\Flowsheet
python flowsheet_tests.py
if %ERRORLEVEL% neq 0 (
    echo [FAIL] Flowsheet tests failed
    set /a FAILED_TESTS+=1
) else (
    echo [OK] Flowsheet tests passed
)
cd ..\..

echo.
echo Testing package installation...
python -c "import pkg_resources; pkg = pkg_resources.get_distribution('PharmaPy'); print('[OK] Package installed:', pkg.project_name, 'v' + pkg.version)"
if %ERRORLEVEL% neq 0 (
    echo [FAIL] Package installation check failed
    set /a FAILED_TESTS+=1
)

echo.
echo ================================================
echo Test Summary
echo ================================================
if %FAILED_TESTS% equ 0 (
    echo  All tests passed!
    exit /b 0
) else (
    echo [FAIL] %FAILED_TESTS% test(s) failed
    exit /b 1
)
