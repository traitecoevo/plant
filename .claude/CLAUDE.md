# CLAUDE.md

This project's architecture and contributor guidance for AI agents lives in
[agents.md](../agents.md). **Read [agents.md](../agents.md) before making changes** —
especially anything touching the C++/R boundary, the RcppR6 interface, or
adding a new model.

Quick pointers:

- Two-layer package: C++ core (`inst/include/plant/`, `src/`) bridged to R
  (`R/`) via RcppR6 + Rcpp. Interface contract: `inst/RcppR6_classes.yml`.
- **Never hand-edit generated files** (`R/RcppR6.R`, `R/RcppExports.R`,
  `src/RcppR6.cpp`, `src/RcppExports.cpp`, `NAMESPACE`, `man/`).
- After interface changes run `make rebuild`; after C++-only changes
  `make compile`; for R-only changes use `devtools::load_all()`.
- Run tests with `make test` or `devtools::test()`.

See [agents.md](../agents.md) for the full architecture, build workflow, and the
process for adding a new strategy/environment model.

## Issue & project-board conventions

Development across `plant`, `plant.assembly`, and `overstorey` is tracked on a
shared [project board](https://github.com/orgs/traitecoevo/projects/5). New issues
are auto-added to the board with status **Backlog** by a workflow, so you do not
need to set status manually.

When opening an issue (including whenever the user asks you to create one), always:

- **Set exactly one type label.** Only three labels exist in these repos — do not
  invent new ones:
  - `bug` — an existing feature not functioning as intended
  - `task` — a discrete task needed for a feature (the default for normal work)
  - `epic` — a new feature or capability, usually an umbrella over several tasks
- **Prefix the title with a theme tag** in square brackets so the board sorts
  cleanly. Reuse an existing theme where it fits; only fall back to `[other]` when
  nothing applies:

  | Tag | Scope |
  |---|---|
  | `[TF24 hydraulics]` | Hydraulics component of the TF24 strategy |
  | `[TF24 allometry]` | Flexible allometry for the TF24 model |
  | `[TF24 nsc]` | Non-structural carbohydrate storage in TF24 |
  | `[acclimation]` | Acclimation of leaf and other traits |
  | `[simplify interface]` | Consistent interface to the plant & plant.assembly models |
  | `[evol assembly]` | Evolutionary assembly linking plant to plant.assembly |
  | `[Env drivers]` | Driving the model with environmental drivers |
  | `[speed]` | Performance — making the model run faster |
  | `[patch variations]` | Multiple patch setups (multi-patch, stochastic metapopulation, continuous patch) |
  | `[AutoDiff]` | Enabling automatic differentiation in plant (e.g. for gradient-based calibration) |
  | `[forecasting]` | Enabling forecasting with the plant model |
  | `[documentation]` | Documenting model capabilities (any of the three repos) |
  | `[other]` | Anything not covered above |

  A title may carry more than one tag when it genuinely spans themes
  (e.g. `[speed] [TF24 hydraulics] …`).

Create issues with `gh issue create -R traitecoevo/plant --title "[tag] …"
--label task` (swap in `bug`/`epic` as appropriate).

## Cross-package context

This repo is part of the **plant family** in the `traitecoevo` org (the core individual-based forest model). For
cross-package orientation — how the family fits together, dependency direction,
source-of-truth rules, and the shared label/board conventions — see
**[`plant-meta`](https://github.com/traitecoevo/plant-meta)** (start with its
`AGENTS.md`). Don't restate family-wide concerns here; link to plant-meta and
keep this file about `plant`-local matters.
