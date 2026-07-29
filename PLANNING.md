# PharmaPy Planning and Release Model

> **Status:** Repository policy when present on the default branch.
>
> This document is the single, version-controlled source for how the
> `PharmaPy-org/PharmaPy` repository plans work and cuts releases. It exists so
> the team can (1) agree on an operating model, (2) generate the GitHub
> milestone, Project views, and release notes *from* a reviewed source rather
> than ad hoc, and (3) verify our own process against an explicit checklist.
>
> Its authority is limited to planning, prioritization, and releases. Coding
> and verification rules belong in `AGENTS.md` if that file is adopted (see
> [§9](#9-adoption-how-these-sections-graduate)). Revisions to either policy
> document go through normal pull-request review.

## Table of contents

1. [Purpose and how to use this document](#1-purpose-and-how-to-use-this-document)
2. [Verified snapshot (2026-07-29)](#2-verified-snapshot-2026-07-29)
3. [Operating model](#3-operating-model)
4. [Release-risk policy](#4-release-risk-policy)
5. [First milestone: `0.1.0a1`](#5-first-milestone-010a1)
6. [Project (#1) changes](#6-project-1-changes)
7. [Decision log](#7-decision-log)
8. [Open questions for maintainers](#8-open-questions-for-maintainers)
9. [Adoption: how these sections graduate](#9-adoption-how-these-sections-graduate)

---

## 1. Purpose and how to use this document

This is a *planning* document, not user-facing docs. Three concrete uses:

- **Verification.** §4–§5 define what "ready to release" means. Before we tag,
  we walk the exit-criteria checklist and it must pass.
- **Generation.** §5 is written to be pasted into a GitHub **milestone**
  description; §6 lists the exact Project view filters to create; §4 is the
  rule we encode into release notes.
- **Discussion.** The team ratifies and revises decisions in §7 by lazy
  consensus (propose a default, set an objection window, silence = assent).

**Guiding principle:** this is a small, asynchronous, academic team whose real
delivery already flows through continuous PR merges. The model below formalizes
*that*, and deliberately avoids ceremony (sprints, mandatory estimates,
blocking committee gates) that the team will not sustain.

## 2. Verified snapshot (2026-07-29)

Queried read-only through the GitHub REST and GraphQL APIs at
**2026-07-29T13:54:25Z**. **Time-sensitive** — re-verify before acting on any
count. Re-running the Project claims requires organization access and a token
with the `read:project` scope.

- 86 open issues; **no milestones, no releases, no Git tags**; `pyproject.toml`
  declares `version = "0.0.1"`.
- Labels on open issues: 62 `correctness`, 60 `status:verified`,
  1 `severity:critical`, 34 `severity:high`, 20 `severity:medium`,
  3 `severity:low`. Severity labels cover 58 of the 62 correctness issues
  (4 correctness issues carry no severity), and every severity-labelled issue
  is a correctness issue.
- The single `severity:critical` is **#68** (adiabatic crystallizer
  energy-balance crash).
- Org **Project #1 "PharmaPy Development"** was created 2026-07-21 and contains
  85 items. Its
  `Priority` values are a mechanical copy of `Severity` (34 High / 20 Medium /
  3 Low / 1 Urgent — identical distribution), the `Size` field is unused
  (0 of 85 items populated), and only #23, #26, #134 sit in any iteration.
- CI (`.github/workflows/ci.yml`) runs core tests on **Python 3.11 only**, while
  `requires-python = ">=3.9"`. Locked pixi installs gate Linux and Windows; the
  Assimulo integration job is `continue-on-error` (informational).
- The codebase dates to 2021 and is published (DOI
  `10.1016/j.compchemeng.2021.107408`); `0.0.1` is a placeholder, not release
  history.

## 3. Operating model

One Project, milestone-anchored, continuous flow. Each GitHub concept has
exactly one meaning:

| Concept | Meaning in this repo |
| --- | --- |
| **Project #1** | The single source of truth: full backlog, active board, and roadmap. We do **not** create a Project per release or subsystem. |
| **Status** | Workflow state only: `Todo` → `In progress` → `Done`. Delivery is controlled by a **work-in-progress (WIP) limit** on `In progress`, not by sprints. Set the cap when a milestone opens using **one active item per active contributor plus one shared review/unblocking slot**; the formula keeps the limit tied to actual team capacity instead of a magic number. |
| **Milestone** | A repository release (or a concrete, externally meaningful outcome). **Exactly one open at a time** because the Project roadmap already holds longer-horizon work; a second open milestone would create a competing release gate. |
| **Epic + sub-issues** | Decomposition of a large initiative, using native issue-type `Epic` and parent/sub-issue links. Current epics: #3, #17, #67, #118. This is the roadmap's backbone. |
| **Priority** | Maintainer delivery order — **re-triaged independently of severity** (see §7 D4). Until re-triaged it is unreliable and is kept out of delivery views. |
| **Severity** (label) | Technical impact of a *defect* only. Does not imply delivery order. |
| **Area** (label) | Affected subsystem. Useful for filtering; not a planning input. |
| **Assignee** | The one person accountable for the item. No duplicate "owner" field. |
| **Target date** | Only for a genuinely externally-dated deliverable. Not a general field. |
| **Iteration** | **Not used as a commitment unit.** (See §7 D3.) Removing it also removes the three-way conflict between iteration dates, `Target date`, and the milestone due date. |
| **Size** | **Not required.** Optional, and only populated for current-milestone items if the team finds it useful. (See §7 D5.) |

Rationale for dropping sprints: the iteration apparatus was created on
2026-07-21 (the Project creation date in §2), its first "iteration" was backdated
and held zero items, and the current board contains work outside the release
while omitting the release's own critical blocker. The team's actual, working
coordination mechanism is per-PR "whoever lands second, rebase" notes —
continuous flow. We formalize that.

## 4. Release-risk policy

The correctness backlog was 62 issues at the §2 snapshot. Blocking a release on
all of it is neither achievable nor necessary. The rule:

- **`severity:critical` → hard blocker.** No release ships with an open
  critical. Each critical must be closed **with a regression test**. *(Today:
  only #68, fixed by PR #106, which adds regression coverage for its four
  sub-defects.)*
- **`severity:high` → not individually blocking, with a promotion rule.** A high
  is *promoted* to a blocker if it either (a) produces a **silent wrong
  numerical result in a documented example/tutorial**, or (b) affects a module
  exercised by the release **smoke/validation suite**. All non-promoted highs
  are **disclosed by number in the release notes**.
- **The gate is behavioral, not lexical.** "No open `severity:critical`" is
  necessary but insufficient, because severity labels do not cover every
  correctness issue. Release also requires the smoke/numerical suite green
  (§5 exit criteria).
- **A named release manager decides promotions**, recorded on the milestone.

## 5. First milestone: `0.1.0a1`

> The block below is written to be pasted into the GitHub milestone
> description once §7 D1/D2 are ratified.

**Milestone title:** `Org transition: installable, documented, CI-verified`
**Release version (Git tag / `pyproject`):** `0.1.0a1` *(alpha — proposed, see §7 D1)*

### Version terminology

This repository uses [PEP 440](https://peps.python.org/pep-0440/) terminology.
The numeric suffix on a pre-release counts successive builds within that phase.

| Term | Meaning in this plan |
| --- | --- |
| **Release** | A uniquely versioned snapshot of the project, represented by an immutable Git tag and published distribution. A release can be a pre-release or a final release. |
| **`0.1.0a1`** | The first **alpha** for the planned `0.1.0` series. It is suitable for early testing while known correctness gaps and planned release work remain. |
| **`0.1.0rc1`** | The first **release candidate**. Planned scope is complete; only fixes for release-blocking findings should separate it from the final release. |
| **`0.1.0`** | The **final release** for this series: the same release segment without a pre-release suffix. “Final” records progression through the release process; it does not by itself claim production validation. |

The proposed progression is `0.1.0a1` → later alphas if findings require them →
`0.1.0rc1` → later release candidates if blockers remain → `0.1.0`. Each later
pre-release increments the numeric suffix for its phase.

**Goal.** Produce the first `PharmaPy-org`-maintained release: the package
installs and its core modules import in a clean, supported environment; core CI
and documentation are reproducible and public; the one critical defect is fixed;
release notes state the remaining known correctness limitations. This is
explicitly an **alpha** — it does not claim numerical production-validation.

### Scope (in) — mapped to tracked issues

Implementation status is a dated snapshot from **2026-07-29T13:54:25Z**; the
issues remain the durable scope records and must be re-queried before copying
this table into a milestone.

| Issue | Implementation status at snapshot | Role |
| --- | --- | --- |
| #68 | PR **#106** open | The critical crystallizer crash + enthalpy basis (hard gate). |
| #134 | PR **#135** open | Solver-free model imports (lazy Assimulo). |
| #130 | PR **#131** open | Reproducible, public documentation (#131 is one step and intentionally does not close #130). |
| #8 | PR **#136** merged | Packaging/metadata modernization, `pharmapy-sim` distribution, pixi environments, install guide, and two-platform locked install matrix. Remaining #8 scope stays on the issue; distribution-name reservation is tracked by #146. |

### Scope (out)

- #7 in its broad form (split out a release-specific smoke/validation issue —
  see §8 Q1).
- #10, the full solver-abstraction initiative.
- Every correctness issue counted in §2 as a single release gate.
- Epic #118 and the complete StateLayout migration.
- Speculative bioreactor, crystallizer-extension, and scheduling initiatives.
- #23 (PR #114) and #26 (PR #115): **in flight, not release scope** unless the
  owner confirms they are intended to ship here (§7 D6). "Work has started" is
  not by itself a reason to include them.

### Entry criteria

- [ ] Every in-scope issue has an owner.
- [ ] #68 has a reproduction and a failing-first regression test.
- [ ] Version scheme (§7 D1) and distribution name (§7 D2) ratified.

### Exit criteria

- [ ] Supported Python versions are explicit and **CI tests exactly what is
      claimed** (either expand CI beyond 3.11 or narrow `requires-python` to
      `>=3.11`). Dependency bounds and environment policy remain sourced in
      [`DEPENDENCIES.md`](DEPENDENCIES.md).
- [x] A clean environment installs the package through the documented pip and
      locked pixi paths. PR #136 added the Linux/Windows smoke-install jobs and
      [`INSTALLATION.md`](INSTALLATION.md); commands and lane semantics are
      sourced in [`TESTING.md`](TESTING.md).
- [ ] Core modules import **without Assimulo**; a solver path without Assimulo
      fails with a clear, localized message (asserted by a test).
- [ ] Core CI and documentation CI are green; the public docs site and
      pull-request previews are verified (#130).
- [ ] **#68 is closed with a regression test**, and the smoke/validation suite
      is green.
- [ ] No open defect is labelled `severity:critical`.
- [ ] Release notes enumerate the open `severity:high` correctness issues by
      number and make no production-validation claim.

### Due-date policy

No due date until the in-scope issues have owners and an agreed capacity
forecast. Once set, the milestone date is a forecast, not a substitute for
issue-level `Target date`.

### Milestone hygiene

Put planned **issues** in the milestone. Do **not** add both an issue and its
linked PR (double-counts progress). Add a PR to the milestone only when the PR
is itself the tracked deliverable with no corresponding issue.

Assignments also live on the milestone **issues**, not in this file or on their
linked PRs. Each in-scope issue has one accountable GitHub assignee. The release
manager (§7 D7) is responsible for assigning unowned work before the milestone
starts and for ensuring the assignee and Project status stay current. When work
is handed off, the current assignee or release manager updates the issue before
the handoff; PR assignees and reviewers describe implementation and review
roles, not milestone ownership.

## 6. Project (#1) changes

**Views** (create/keep exactly these; remove the rest):

| View | Layout | Filter |
| --- | --- | --- |
| Backlog | Table | `is:open` |
| Board | Board (by Status) | `is:open` — enforce the WIP limit on `In progress` |
| Release: 0.1.0a1 | Table | `milestone:"Org transition: installable, documented, CI-verified"` |
| Correctness defects | Table | `label:correctness` sorted by severity |
| Unassigned high-priority | Table | `no:assignee (label:severity:critical OR label:severity:high)` |
| Roadmap | Roadmap | `type:Epic` (not every Project item) |

Remove the `Current iteration` view (no iterations).

**Fields:** retire `Iteration` and `Size` from the model; keep Status, Priority
(re-triaged), Severity/Area (labels), Milestone, Assignee, and a sparingly-used
`Target date`.

**Automations:** keep auto-add of open issues and archive-on-close; add
auto-set `Status: Done` when the closing PR merges. Stop mirroring Severity into
Priority.

Operational follow-ups are recorded on their owning issues rather than as
one-shot checkboxes here: #68 owns hard-blocker status, while #134 owns its
assignee and Project status.

## 7. Decision log

Ratify by lazy consensus. Each ratification announcement names its proposed
default and exact start/end timestamps for a **one-week objection window**;
silence = assent. One week is used because it spans a full workweek for
asynchronous contributors. The announcement belongs in the relevant pull
request or tracking issue, not as a transient deadline in this file. The owner
makes the call if consensus is unclear.

| # | Decision | Proposed default | Reversible? | Owner |
| --- | --- | --- | --- | --- |
| **D1** | First release version | `0.1.0a1` (alpha) | Hard once tagged/DOI'd | maintainer |
| **D2** | Distribution name on the index | `pharmapy-sim` (recheck at publish, #146) | Hard once published | maintainer |
| **D3** | Timeboxed iterations | **Drop**; use WIP-limited continuous flow | Yes | team |
| **D4** | Priority model | Re-triage independently of severity (or collapse to `Now/Next/Later`) | Yes | maintainer |
| **D5** | `Size` field | Not required; optional for near-term items only | Yes | team |
| **D6** | #23/#26 in this release? | No unless the owner confirms they ship here | Yes | issue owner |
| **D7** | Release manager | The maintainer who creates the milestone acts as release manager unless the milestone names another volunteer | Yes | maintainer |

**Version comparison for D1:**

- `0.0.2` — too timid; perpetuates the misleading `0.0.1` lineage and understates
  a deliberate first org release.
- `0.1.0` — honest about the capability step, but reads as a stable feature
  baseline, which 34 open high-severity correctness defects contradict.
- **`0.1.0a1`** *(proposed)* — signals "this is the target shape; correctness
  work continues," and allows the alpha → release-candidate → final progression
  defined in §5 as the backlog burns down.
- A non-versioned milestone name is used regardless (above); the *tag* still
  needs a number, so pair it with `0.1.0a1`.

## 8. Open questions for maintainers

These need a human decision; do not invent answers.

1. **Smoke/validation suite.** [`TESTING.md`](TESTING.md) already defines the
   core pytest lane, locked pixi install matrix, and informational Assimulo
   lane. Which of those existing lanes gate the release, and what numerical
   assertions must be added on top? Should #7 be split into a release-specific
   testing issue?
2. **Supported matrix.** Honor `>=3.9` (and test it in CI) or narrow to
   `>=3.11`?
3. **Citable release.** Given the 2021 paper and DOI, do we want a
   `CITATION.cff` + Zenodo archive for the eventual `0.1.0`? (Not required for
   the alpha.)
4. **Safety/validation bar.** Do any users consume PharmaPy numerical output for
   real process decisions? If so, the acceptable correctness bar may exceed "no
   open critical," and the release-notes disclaimer must reflect it.
5. **Release cadence after 0.1.0.** Patch milestones (`0.1.x`) for bounded
   correctness fixes; a later minor for meaningful API/architecture change. No
   future milestones created until their scope is credible.

## 9. Adoption: how these sections graduate

On merge, this document becomes the reference. As decisions are ratified:

- §4 (release-risk policy) and §5 due-date/hygiene → the durable content of a
  future `RELEASING.md`.
- §3 (operating model) → the planning half of a future `CONTRIBUTING.md` / the
  Project README.
- §5 (milestone block) → pasted into the created GitHub milestone.
- §6 → applied to Project #1.

The authority split with PR #132 is deliberate: `PLANNING.md` owns work
planning, prioritization, and releases; the proposed `AGENTS.md` owns coding and
verification rules. A future contributor-facing `CONTRIBUTING.md` should
summarize and point to both rather than duplicate either.

### Maintenance responsibility and refresh triggers

Any maintainer may propose policy changes through a pull request. The named
release manager owns milestone-specific maintenance: keeping the milestone
description synchronized with §5 and confirming the release gates before a tag.

Maintenance is event-driven rather than scheduled on a fixed calendar. Refresh
the dated §2 snapshot and §5 implementation status together:

- before creating or rolling over a milestone;
- before tagging any pre-release or final release; and
- when a ratified decision changes scope, the operating model, or a release
  gate.

Changing an individual work assignment does not require a `PLANNING.md` edit:
the issue assignee and Project status are the live operational record. A policy
change to who owns assignments does require review here. Do not advance the
snapshot timestamp unless every snapshot claim is re-queried in the same pass.

This file stays as the living planning record; the generated artifacts
(milestone, views, `RELEASING.md`) are downstream of it.
