# Dependency Policy

The default install keeps only the core runtime dependencies. Assimulo is an
optional solver dependency because it is difficult to build on common pip-only
environments and is only needed for the solver-backed unit-operation models.

Current bounds:

`pyproject.toml` is the source of truth for install metadata. The
`requirements.txt` and `requirements-assimulo.txt` files are convenience
mirrors for pip-based workflows and should stay synchronized with
`[project.dependencies]` and the `assimulo` optional extra.

- `numpy>=1.22`: the lower bound gives a modern baseline for supported Python
  versions while allowing the NumPy 2 line to be exercised by CI and Dependabot.
- `scipy>=1.9`: the previous SciPy upper cap was tied to
  `scipy.integrate.simps`, which was removed in SciPy 1.14. PharmaPy now uses
  `scipy.integrate.simpson`, so current SciPy releases can be tested directly.
- `matplotlib>=3.5` and `pandas>=1.5`: lower bounds avoid unbounded old
  releases while leaving current releases open for Dependabot and CI.
- `cython>=0.29,<3` and `assimulo==3.4.3`: these are confined to the
  `assimulo` optional extra until the Assimulo integration job proves a broader
  solver stack works. Dependabot intentionally ignores Assimulo updates and
  Cython major-version updates until that broader solver stack is verified.
