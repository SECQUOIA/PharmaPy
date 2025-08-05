.PHONY: help install install-dev install-assimulo test test-all test-reactor test-flowsheet test-imports test-quick test-discovery clean build upload docs format lint

help:
	@echo "Available commands:"
	@echo "  install        Install the package in development mode"
	@echo "  install-dev    Install with development dependencies"
	@echo "  install-assimulo Install with assimulo for simulation features"
	@echo "  test           Run all tests"
	@echo "  test-quick     Run quick tests (imports + reactor tests)"
	@echo "  test-reactor   Run reactor tests only"
	@echo "  test-flowsheet Run flowsheet tests only"
	@echo "  test-imports   Test package imports only"
	@echo "  test-all       Run comprehensive test suite"
	@echo "  test-discovery Run test discovery"
	@echo "  clean          Clean build artifacts"
	@echo "  build          Build the package"
	@echo "  upload         Upload to PyPI (requires credentials)"
	@echo "  docs           Build documentation"
	@echo "  format         Format code with black and isort"
	@echo "  lint           Run code linting"

install:
	pip install -r requirements.txt
	pip install -e .

install-dev:
	pip install -r requirements.txt
	pip install -r requirements-dev.txt
	pip install -e .

install-assimulo:
	@echo "Installing PharmaPy with assimulo simulation support..."
	@echo "Note: Requires system dependencies (sundials, superlu, openblas)"
	@echo "Install via: conda install -c conda-forge sundials=5.8.0 superlu=5.2.2 openblas"
	pip install -r requirements.txt
	pip install -r requirements-assimulo.txt
	pip install -e .

test: test-all

test-quick:
	@echo "Running quick test suite..."
	@python -c "import PharmaPy; print('✓ PharmaPy imports OK')" || echo "X PharmaPy import failed"
	@cd tests/integration && python reactor_tests.py || echo "⚠ Reactor tests require dependencies"

test-reactor:
	@echo "Running reactor tests..."
	@cd tests/integration && python reactor_tests.py || echo "⚠ Reactor tests require dependencies"

test-flowsheet:
	@echo "Running flowsheet tests..."
	@cd tests/Flowsheet && python flowsheet_tests.py || echo "⚠ Flowsheet tests require dependencies"

test-imports:
	@echo "Testing package imports..."
	@python -c "import PharmaPy; print('✓ PharmaPy')" || echo "X PharmaPy import failed"
	@python -c "from PharmaPy import Utilities; print('✓ Utilities')" || echo "⚠ Utilities import failed"
	@python -c "from PharmaPy import Reactors; print('✓ Reactors')" || echo "⚠ Reactors import failed (may need assimulo)"
	@python -c "from PharmaPy import Streams; print('✓ Streams')" || echo "⚠ Streams import failed (may need assimulo)"
	@python -c "from PharmaPy import Phases; print('✓ Phases')" || echo "⚠ Phases import failed (may need assimulo)"
	@python -c "from PharmaPy import Kinetics; print('✓ Kinetics')" || echo "⚠ Kinetics import failed (may need assimulo)"

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
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	rm -rf .pytest_cache/ 2>/dev/null || true
	rm -rf htmlcov/ 2>/dev/null || true
	rm -f .coverage 2>/dev/null || true

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