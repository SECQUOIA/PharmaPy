.PHONY: help install install-dev test clean build upload docs format lint

help:
	@echo "Available commands:"
	@echo "  install      Install the package in development mode"
	@echo "  install-dev  Install with development dependencies"
	@echo "  test         Run tests"
	@echo "  clean        Clean build artifacts"
	@echo "  build        Build the package"
	@echo "  upload       Upload to PyPI (requires credentials)"
	@echo "  docs         Build documentation"
	@echo "  format       Format code with black and isort"
	@echo "  lint         Run code linting"

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

test:
	cd tests && python reactor_tests.py
	@echo "Tests completed successfully!"

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	find . -type d -name __pycache__ -delete
	find . -type f -name "*.pyc" -delete

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
	cd tests && python reactor_tests.py

clean-win:
	if exist build rmdir /s /q build
	if exist dist rmdir /s /q dist
	if exist *.egg-info rmdir /s /q *.egg-info
	for /d /r . %%d in (__pycache__) do @if exist "%%d" rmdir /s /q "%%d"
	del /s /q *.pyc 2>nul
