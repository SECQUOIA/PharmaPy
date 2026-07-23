# Installing PharmaPy

PharmaPy has a pip-installable core and an optional Assimulo solver backend.
Choose the installation path that matches the work you need to do:

- Use **pip** for the core library and solver-free development.
- Use **pixi** for a reproducible development environment or any work that
  requires Assimulo.

pixi is a development and solver-environment tool. It is not a PharmaPy runtime
dependency and does not replace the package metadata used by pip.

The installation commands intentionally produce different environments:

| Command | Assimulo included? | Intended use |
| --- | --- | --- |
| `python -m pip install .` | No | Core library for users |
| `pixi install --locked` | No | Reproducible core development |
| `pixi install -e assimulo --locked` | Yes | Solver-backed development and use |
| `pixi install --all --locked` | Yes, in the separate `assimulo` environment | Maintainer verification of both environments |

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

## Why pixi is used

pip remains the primary installation mechanism for the PharmaPy package. pixi
is used for the development and solver environment because Assimulo wraps
compiled SUNDIALS, Cython, C, and Fortran components. The
[Assimulo release on PyPI](https://pypi.org/project/Assimulo/) does not provide
a compatible wheel for PharmaPy's supported Python versions, whereas
[conda-forge provides compiled Assimulo packages](https://anaconda.org/conda-forge/assimulo).

pixi can resolve those conda-forge binaries, install the local PharmaPy package
editably, and record the complete result in `pixi.lock`. Keeping pixi
configuration under `[tool.pixi]` means that it does not become a dependency of
the PharmaPy wheel and does not change how the package is built or published to
PyPI.

Plain `environment.yml` remains a simpler fallback, but it does not pin the
complete transitive environment. The lockfile is the reason to retain pixi; if
the project stops testing and maintaining that lockfile, pixi would add
maintenance cost without providing its intended reproducibility benefit.

## Reproducible installation with pixi

The committed `pixi.lock` currently supports `linux-64`. Install pixi using its
[official installation instructions](https://pixi.prefix.dev/latest/installation/),
then clone and enter the repository:

```bash
git clone https://github.com/PharmaPy-org/PharmaPy.git
cd PharmaPy
```

The available environments are:

- `default`: core dependencies without Assimulo
- `assimulo`: the full solver stack with Assimulo 3.4.3 and SUNDIALS

For work that requires Assimulo, install its locked environment explicitly:

```bash
pixi install -e assimulo --locked
```

This command installs both PharmaPy and Assimulo in the named environment.
Verify that both imports work:

```bash
pixi run -e assimulo python -c "import PharmaPy, assimulo"
```

For core-only development, install the default environment:

```bash
pixi install --locked
```

Maintainers can install both environments together:

```bash
pixi install --all --locked
```

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

## What is tested

The following installation paths were run successfully on `linux-64` when this
guide was introduced:

- A fresh Python 3.11 virtual environment ran `python -m pip install .`,
  `python -m pip check`, and `import PharmaPy`. Wheel and source-distribution
  builds also passed.
- `pixi install --all --locked` installed both locked environments.
- The core pixi test lane passed.
- The Assimulo pixi test lane passed and imported both `PharmaPy` and Assimulo
  3.4.3.

Continuous integration currently exercises the core package through an
editable pip install and exercises Assimulo through `environment.yml`. It does
not yet execute the exact pixi installation command. Until a pixi CI lane is
added, maintainers changing dependencies should run:

```bash
pixi lock --check
pixi install --all --locked
pixi run test
pixi run -e assimulo test-assimulo
```

Keeping the pip-based core CI lane alongside a locked pixi Assimulo lane proves
both promises: PharmaPy remains pip-installable, and the compiled solver
environment remains reproducible.

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
