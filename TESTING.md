# Testing

Run the core pytest slice without the optional Assimulo solver stack:

```bash
python -m pip install -e ".[test]"
python -m pytest --collect-only
python -m pytest tests/ -m "not assimulo"
```

Run the Assimulo-backed integration tests with conda-forge:

```bash
conda env create -f environment.yml
conda activate pharmapy-assimulo
python -m pip install -e . --no-deps
python -m pytest tests/ -v -m assimulo
```

`python -m pip install -e ".[assimulo]"` installs only pip-available build
helpers for source workflows; it does not install Assimulo itself. Use the
conda environment above for the supported Assimulo test path.

The GitHub Actions Assimulo job is intentionally informational at first. It
uses `continue-on-error: true` so failures in the external solver stack do not
block the core test job while the environment is stabilized.
