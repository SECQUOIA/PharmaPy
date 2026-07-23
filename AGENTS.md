# PharmaPy Agent Guidelines

This file is the canonical project guidance for coding agents. Read it before
making changes. `CLAUDE.md` and `.github/copilot-instructions.md` load or point
to this file so that Codex, Claude Code, and GitHub Copilot share one policy.

## Project vision

PharmaPy is scientific software for modeling and simulating pharmaceutical
manufacturing processes. Changes must preserve the physical meaning of the
models and make simulations trustworthy, reproducible, and understandable to
chemical-engineering users.

Prioritize, in order:

1. Physical and numerical correctness.
2. Clear units, bases, assumptions, and APIs.
3. Backward compatibility and reproducibility.
4. Maintainability.
5. Performance, when supported by measurements.

## Working principles

- Read the relevant implementation, tests, and documentation before editing.
- Make the smallest coherent change that solves the stated problem. Do not
  combine unrelated refactors with a fix or feature.
- Preserve public APIs and established model behavior unless the task explicitly
  requires a breaking change. Document any intentional compatibility impact.
- Never hide a physical-basis conversion or silently change a model assumption.
- Treat test data and expected numerical results as scientific evidence; do not
  update them merely to make a failing test pass.
- Do not commit generated files, local environments, credentials, or machine-
  specific configuration.

## Python documentation

- Every new or modified function and method, including private helpers, must
  have a NumPy-style docstring. When touching an undocumented function, bring
  its docstring into compliance as part of the change.
- Include `Parameters`, `Returns` or `Yields`, and `Raises` sections when they
  apply. Include `Notes`, `Warnings`, or `Examples` when assumptions, equations,
  side effects, or usage would otherwise be unclear.
- Document the units, physical basis, expected shape, and meaning of every
  scientific input and output. For example: `temperature : float` followed by
  `Temperature of the liquid phase [K].`
- Class docstrings must describe the model's purpose and important assumptions;
  method docstrings must document method-specific behavior rather than relying
  on the class description.
- Comments should explain intent, assumptions, or non-obvious mathematics. Do
  not narrate code that is already clear.

## Units and physical quantities

- Every new or modified variable that stores a numeric physical or model
  quantity must have an adjacent unit comment in square brackets at its first
  definition. Use SI-derived notation where practical:

  ```python
  particle_diameter = 2.5e-4  # [m]
  mass_flow = inlet.mass_flow  # [kg/s]
  conversion = reacted_moles / initial_moles  # [-]
  ```

- Use `[-]` for dimensionless quantities. Pure bookkeeping values such as
  indices, booleans, names, and object references do not require a unit comment.
- For arrays, mappings, or tuples with heterogeneous units, state the unit of
  each component in the adjacent comment or NumPy-style docstring.
- Function parameters and return values must state units in their docstrings.
  Attributes, constants, intermediate calculations, and test fixtures must use
  adjacent unit comments.
- Repeat or update the unit comment whenever a conversion changes the unit or
  basis. Use explicit variable names for conversions instead of reusing a name
  with a different meaning.
- State the basis as well as the unit when ambiguity is possible, such as
  `[J/mol]`, `[J/kg]`, `[mol/s]`, or `[kg/s]`. Never mix mass, molar, volume, or
  particle-number bases implicitly.
- Check dimensional consistency in every equation. Centralize conversions and
  avoid unexplained numeric conversion factors or magic constants.

## Numerical and API conventions

- Use descriptive, domain-relevant names. Avoid one-letter names except for
  conventional short-lived mathematical indices or coordinates.
- Add type hints to new or modified public APIs and to internal code when they
  clarify array shapes, optional values, or return structure.
- Preserve scalar-versus-array behavior and document expected NumPy shapes.
- Validate inputs at public boundaries and raise specific exceptions with
  actionable messages for invalid units, shapes, ranges, or model state.
- Use explicit, justified tolerances for floating-point comparisons. Tests
  should use `pytest.approx` or NumPy testing helpers rather than exact equality
  for calculated floating-point values.
- Handle zero, negative, empty, NaN, and infinite values deliberately when they
  are plausible at a model boundary.
- Prefer existing project abstractions and NumPy/SciPy operations over duplicate
  helpers. Add a dependency only when the benefit justifies the maintenance cost.

## Testing and verification

- Add or update tests for every behavior change and bug fix. A regression test
  should fail without the fix and pass with it.
- Include nominal behavior plus relevant boundary, dimensional, and failure
  cases. Assert units or physical basis indirectly through known balances or
  conversions when practical.
- Use the markers defined in `pytest.ini`: `unit`, `integration`, `slow`, and
  `assimulo`. Do not make core tests depend on the optional Assimulo stack.
- Run the narrowest relevant tests while developing, then run the core suite
  before handoff:

  ```bash
  python -m pytest tests/ -m "not assimulo"
  ```

- For solver-backed changes, also use the conda environment and Assimulo test
  command documented in `TESTING.md`. If that environment is unavailable, say
  exactly which verification was not run and why.
- Keep `pyproject.toml`, dependency mirrors, and `DEPENDENCIES.md` synchronized
  when changing dependencies.

## Documentation and handoff

- Update user documentation and examples when public behavior, model equations,
  inputs, outputs, units, or installation steps change.
- In the final change summary, explain the scientific or behavioral effect,
  identify tests run, and call out assumptions, compatibility concerns, and
  unverified optional-solver paths.
- If a request conflicts with these guidelines, surface the conflict and obtain
  maintainer direction instead of silently ignoring the policy.
