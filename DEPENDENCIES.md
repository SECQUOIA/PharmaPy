# Dependency Policy

The default install keeps only the core runtime dependencies. Assimulo is an
optional solver dependency because it is difficult to build on common pip-only
environments, is not published on PyPI above version 3.0, and is only needed for
the solver-backed unit-operation models. For solver-backed runs, use the pixi
environments defined in `pyproject.toml` (`pixi run -e assimulo ...`), the
conda-forge environment in `environment.yml`, or a local source build of
Assimulo.

pixi is the recommended developer/CI environment manager: its `[tool.pixi]`
configuration in `pyproject.toml` provides a `default` (core, Assimulo-free)
and an `assimulo` (full backend) environment backed by a committed `pixi.lock`,
which is the only reproducible way to obtain Assimulo since it is conda-forge
only. `[project]` remains the source of truth for pip users; the pixi
dependency bounds mirror it and must stay synchronized.

Current bounds:

`pyproject.toml` is the source of truth for install metadata. The
`requirements.txt` and `requirements-assimulo.txt` files are convenience
mirrors for pip-based workflows and should stay synchronized with
`[project.dependencies]` and the pip-installable entries in the `assimulo`
optional extra.

- `numpy>=1.22`: the lower bound gives a modern baseline for supported Python
  versions while allowing the NumPy 2 line to be exercised by CI and Dependabot.
- `scipy>=1.9`: the previous SciPy upper cap was tied to
  `scipy.integrate.simps`, which was removed in SciPy 1.14. PharmaPy now uses
  `scipy.integrate.simpson`, so current SciPy releases can be tested directly.
- `matplotlib>=3.5` and `pandas>=1.5`: lower bounds avoid unbounded old
  releases while leaving current releases open for Dependabot and CI.
- `cython>=0.29,<3`: kept in the `assimulo` optional extra as a pip-installable
  build helper for source-install workflows. The extra does not install
  Assimulo itself.
- `assimulo==3.4.3`: the supported Assimulo version for the conda-forge
  environment and source-build workflows. It is documented here instead of
  pinned in pip metadata because that release is not available from PyPI.
  Dependabot intentionally ignores Assimulo updates and Cython major-version
  updates until a broader solver stack is verified.
