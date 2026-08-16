# Contributing to PharmaPy

Thank you for contributing to PharmaPy. This document describes the
development workflow: how changes are proposed, kept in sync, reviewed, and
merged. The engineering standards themselves — physical correctness, units,
docstrings, testing, and verification — live in [AGENTS.md](AGENTS.md) and
are binding for every contribution, whether written by a person or a coding
agent.

## Where the rules live

- [AGENTS.md](AGENTS.md) — canonical coding, units, documentation, and
  testing standards for all changes.
- [TESTING.md](TESTING.md) — test environments and how to run each lane.
- [DEPENDENCIES.md](DEPENDENCIES.md) — dependency policy and environments.
- [INSTALLATION.md](INSTALLATION.md) — installation instructions.

## Development setup

PharmaPy uses [pixi](https://pixi.sh) for reproducible development
environments backed by conda-forge:

```bash
pixi run test                       # core test lane (Assimulo-free)
pixi run -e assimulo test-assimulo  # solver-backed integration lane
pixi run docs                       # Sphinx docs build (warnings are errors)
```

`python -m pytest tests/ -m "not assimulo"` is the equivalent core command in
an already-provisioned environment; see [TESTING.md](TESTING.md).

## Pull request basics

- Branch from `master`. Keep one concern per pull request; do not combine
  refactors with fixes or features.
- Every behavior change carries tests, and CI must be green before review is
  requested (see [AGENTS.md](AGENTS.md) for the testing and verification
  bar).
- Reference the issue a PR addresses; use `Closes #N` only when every
  acceptance criterion of the issue is satisfied.

## Integration contract for parallel workstreams

A *workstream* is any stream of development expected to run for multiple
weeks in parallel with `master` — a new capability module, an architecture
refactor, or a large migration. The clauses below exist because long-lived
branches that drift silently are far more expensive to merge than the weekly
discipline they replace. They apply from the stream's first commit.

1. **Branch from `master` and open a draft pull request within the first
   week** of work, so the stream is visible from day one.
2. **Track every stream on the team project board**
   ([PharmaPy Development](https://github.com/orgs/PharmaPy-org/projects/1))
   with an owner and a current status.
3. **Update from `master` at least weekly**, by rebase or by an append-only
   merge from `master`. A stream must never sit more than two weeks behind
   `master`, and work is sliced so each PR can merge within about two weeks
   of opening.
4. **Tests accompany every PR** in the stream, per [AGENTS.md](AGENTS.md).
   Parallel streams get no exemption: a stream that merges without tests
   transfers its verification debt to whoever integrates after it.
5. **Unit comments on every new physical quantity** at first definition, for
   example `mass_flow = inlet.mass_flow  # [kg/s]`, and NumPy-style
   docstrings stating units, bases, and assumptions, per
   [AGENTS.md](AGENTS.md).
6. **No `_old`, `_refactor`, scratch, or editor-configuration files** in a
   production branch. Superseded code is deleted, not renamed; history is
   preserved by git.
7. **Declare the merge order when two streams touch the same modules.** The
   newer, smaller, better-tested stream reconciles onto the other. When no
   order has been agreed, an architecture refactor merges before capability
   streams that build on the affected modules.
8. **Acknowledge this contract before the first commit.** Each incoming
   workstream posts a short acknowledgement on its tracking issue confirming
   the clauses above.

Maintainers may pause a stream that repeatedly falls out of contract (for
example, two consecutive missed weekly syncs) and land what exists before
the stream continues.
