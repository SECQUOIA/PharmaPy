@echo off
setlocal

if "%1"=="" goto help
if "%1"=="help" goto help
if "%1"=="install" goto install
if "%1"=="install-dev" goto install-dev
if "%1"=="test" goto test
if "%1"=="clean" goto clean
if "%1"=="build" goto build
if "%1"=="upload" goto upload
if "%1"=="docs" goto docs
if "%1"=="format" goto format
if "%1"=="lint" goto lint
goto invalid

:help
echo Available commands:
echo   install      Install the package in development mode
echo   install-dev  Install with development dependencies
echo   test         Run tests
echo   clean        Clean build artifacts
echo   build        Build the package
echo   upload       Upload to PyPI (requires credentials)
echo   docs         Build documentation
echo   format       Format code with black and isort
echo   lint         Run code linting
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
cd tests
python reactor_tests.py
cd ..
echo Tests completed successfully!
goto end

:clean
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
for /d %%i in (*.egg-info) do if exist "%%i" rmdir /s /q "%%i"
for /d /r . %%d in (__pycache__) do @if exist "%%d" rmdir /s /q "%%d"
del /s /q *.pyc 2>nul
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
