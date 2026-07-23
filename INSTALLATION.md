# Installing PharmaPy

PharmaPy has a pip-installable core and an optional Assimulo solver backend.
Choose the installation path that matches the work you need to do:

- Use **pip** for the core library and solver-free development.
- Use **pixi** for a reproducible development environment or any work that
  requires Assimulo.

pixi is a development and solver-environment tool. It is not a PharmaPy runtime
dependency and does not replace the package metadata used by pip.

## Package and import names

The planned PyPI distribution name is **`pharmapy-sim`**. The Python import name
remains **`PharmaPy`**:

```python
import PharmaPy
```

There is not yet a PharmaPy process-simulation release on PyPI. Until one is
published, install this project from a source checkout as described below. Do
not run `pip install pharmapy`: that name belongs to an unrelated project.

## Core installation with pip

The core installation supports Python 3.9 or newer and does not install
Assimulo.

Clone the repository and enter it:

```bash
git clone https://github.com/PharmaPy-org/PharmaPy.git
cd PharmaPy
```

Create and activate a virtual environment:

```bash
python -m venv .venv
```

On Linux or macOS:

```bash
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install PharmaPy:

```bash
python -m pip install --upgrade pip
python -m pip install .
```

Verify the installation:

```bash
python -c "import PharmaPy; print(PharmaPy.__file__)"
```

For editable development and the core test dependencies, use:

```bash
python -m pip install -e ".[test]"
python -m pytest tests/ -m "not assimulo"
```

After the first `pharmapy-sim` release is published, the core installation
command will be:

```bash
python -m pip install pharmapy-sim
```

The `assimulo` optional extra contains only the Cython helper used by source
build workflows; it cannot install the current supported Assimulo release from
PyPI.

## Reproducible installation with pixi

The committed `pixi.lock` currently supports `linux-64`. Install pixi using its
[official installation instructions](https://pixi.prefix.dev/latest/installation/),
then clone and enter the repository:

```bash
git clone https://github.com/PharmaPy-org/PharmaPy.git
cd PharmaPy
```

Install both locked environments without changing the lockfile:

```bash
pixi install --all --locked
```

The environments are:

- `default`: core dependencies without Assimulo
- `assimulo`: the full solver stack with Assimulo 3.4.3 and SUNDIALS

Run the matching test lanes:

```bash
pixi run test
pixi run -e assimulo test-assimulo
```

Run Python or a script in either environment:

```bash
pixi run python -c "import PharmaPy"
pixi run -e assimulo python -c "import PharmaPy, assimulo"
pixi run -e assimulo python path/to/script.py
```

You can also start an activated shell with `pixi shell` or
`pixi shell -e assimulo`.

macOS and Windows are not currently present in the committed lockfile.
Maintainers extending platform support must add the platform under
`[tool.pixi.workspace].platforms`, regenerate `pixi.lock`, and verify both test
lanes on that platform.

## Manual conda fallback for Assimulo

`environment.yml` remains available as an unlocked conda-forge fallback:

```bash
conda env create -f environment.yml
conda activate pharmapy-assimulo
python -m pip install -e . --no-deps
python -m pytest tests/ -m assimulo
```

Prefer pixi when possible because its committed lockfile records the complete
resolved environment.
