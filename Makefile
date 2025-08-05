.PHONY: help install install-dev test test-all test-reactor test-flowsheet test-imports test-quick clean build upload docs format lint

help:
	@echo "Available commands:"
	@echo "  install        Install the package in development mode"
	@echo "  install-dev    Install with development dependencies"
	@echo "  test           Run all tests"
	@echo "  test-quick     Run quick tests (imports + reactor tests)"
	@echo "  test-reactor   Run reactor tests only"
	@echo "  test-flowsheet Run flowsheet tests only"
	@echo "  test-imports   Test package imports only"
	@echo "  test-all       Run comprehensive test suite"
	@echo "  clean          Clean build artifacts"
	@echo "  build          Build the package"
	@echo "  upload         Upload to PyPI (requires credentials)"
	@echo "  docs           Build documentation"
	@echo "  format         Format code with black and isort"
	@echo "  lint           Run code linting"

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

test: test-all

test-quick:
	@echo "Running quick test suite..."
	python -c "import PharmaPy; print('✅ PharmaPy imports OK')"
	@cd tests/integration && python reactor_tests.py

test-reactor:
	@echo "Running reactor tests..."
	@cd tests/integration && python reactor_tests.py

test-flowsheet:
	@echo "Running flowsheet tests..."
	@cd tests/Flowsheet && python flowsheet_tests.py

test-imports:
	@echo "Testing package imports..."
	python -c "import PharmaPy; print('✅ PharmaPy')"
	python -c "from PharmaPy import Reactors; print('✅ Reactors')"
	python -c "from PharmaPy import Streams; print('✅ Streams')" 
	python -c "from PharmaPy import Phases; print('✅ Phases')"
	python -c "from PharmaPy import Kinetics; print('✅ Kinetics')"

test-all:
	@echo "Running comprehensive test suite..."
	python run_tests.py

test-discovery:
	@echo "Running test discovery..."
	python test_discovery.py

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	find . -type d -name __pycache__ -delete
	find . -type f -name "*.pyc" -delete
	rm -rf .pytest_cache/
	rm -rf htmlcov/
	rm -f .coverage

build: clean
	python -m build

upload: build
	twine upload dist/*

docs:
	cd doc && make html

format:
	black PharmaPy/
	isort PharmaPy/

lint:
	flake8 PharmaPy/

# Windows-compatible commands
install-win:
	pip install -e .

install-dev-win:
	pip install -e ".[dev]"

test-win:
	run_tests.bat

test-quick-win:
	python -c "import PharmaPy; print('✅ PharmaPy imports OK')"
	cd tests\integration && python reactor_tests.py

clean-win:
	if exist build rmdir /s /q build
	if exist dist rmdir /s /q dist
	if exist *.egg-info rmdir /s /q *.egg-info
	for /d /r . %%d in (__pycache__) do @if exist "%%d" rmdir /s /q "%%d"
	del /s /q *.pyc 2>nul
	if exist .pytest_cache rmdir /s /q .pytest_cache
	if exist htmlcov rmdir /s /q htmlcov
	if exist .coverage del .coverage
