---
name: plant-new-strategy
description: >-
  Scaffold a new strategy (and optionally environment) model in the plant
  C++/R package, then wire it through the build. Use when asked to add a new
  model/strategy, create a strategy from an existing one (e.g. clone FF16),
  add a variant strategy that reuses an existing environment, or "use the
  scaffolder". Covers running scripts/new_strategy_scaffolder.R, the
  own-environment vs reuse-environment modes (issue #274), implementing the
  biology, and the make rebuild + test loop. Does NOT cover changing the
  numerics of an existing model.
---

# Adding a new strategy/environment model to plant

A **strategy** is plant's biological model for one species (growth, mortality,
reproduction, competition). An **environment** is the resource field a strategy
reads (e.g. the light profile). The C++ core is templated on the
`<Strategy, Environment>` pair, so adding a model means teaching every templated
class and every R dispatch table about the new pair. The scaffolder automates
the boilerplate; **you still write the biology by hand.**

Read [agents.md](../../../agents.md) §7 first. Never hand-edit generated files
(`R/RcppR6.R`, `R/RcppExports.R`, `src/RcppR6.cpp`, `src/RcppExports.cpp`,
`NAMESPACE`, `man/`).

## Two modes

| Mode | When | What it generates |
|---|---|---|
| **Own environment** (default) | The new model needs a different resource field, or you'll change the environment. | A full `<Name>_Strategy` **and** `<Name>_Environment` (new header, yml block, bindings, dispatch). |
| **Reuse environment** (`environment=`) | A variant strategy that competes in an existing environment — e.g. an FF16-variant that shares `FF16_Environment` (issue #274). | Only `<Name>_Strategy`. No `<Name>_Environment` files/yml/bindings are created; the pair reuses the named environment everywhere. |

Naming convention: two initials + year (e.g. `FF16`), optional single-letter
suffix for a minor variant (e.g. `FF16r`).

## Step 1 — scaffold

From the package root:

```r
source("scripts/new_strategy_scaffolder.R")

# Own environment, cloned from FF16:
create_strategy_scaffold("XX24", template_strategy = "FF16")

# Variant reusing an existing environment (no XX24r_Environment generated):
create_strategy_scaffold("FF16r", template_strategy = "FF16", environment = "FF16")
```

`template_strategy` is the model whose files are cloned and renamed (`FF16`,
`K93`, or `TF24`). `environment` (optional) names an existing model whose
`<env>_Environment` will be reused instead of generating a new one.

The scaffolder edits/creates:

- `R/<name>.R`, `src/<name>_strategy.cpp`, `src/<name>_node.cpp`,
  `inst/include/plant/models/<name>_strategy.h`, `tests/testthat/test-strategy-<name>.R`
  (plus `inst/include/plant/models/<name>_environment.h` in own-env mode);
- the `concrete:`/`templates:` pairs and yml blocks in `inst/RcppR6_classes.yml`;
- the `#include` in `inst/include/plant.h`;
- the export in `src/individual_runner.cpp`;
- the dispatch `switch()`es in `R/strategy_support.R` (make_hyperpar,
  param_hyperpar, hyperpar, environment_type, Environment, expand_state,
  node_schedule_default, make_node_schedule);
- the three lists in `tests/testthat/helper-plant.R`.

The scaffolder uses **fixed-string anchors** keyed off the template model's
lines and errors loudly (`anchor not found`) if the codebase has drifted — if
that happens, the anchor strings in `scripts/new_strategy_scaffolder.R` need
updating to match the current dispatch tables.

## Step 2 — implement the biology

Edit the generated strategy (and environment, if own) to define the model:
growth/mortality/reproduction rates, `compute_rates()`, the competition→
environment mapping, and `state_names()`/`state_size()`. Use the cloned template
as the reference. The full Kohyama-1993 walkthrough — the C++ rate functions,
environment mapping, `compute_rates`, header declarations, yml block and
hyperpar — is in [worked-example-k93.md](worked-example-k93.md). Update the
`<Name>_Strategy` yml `list:` block to match the parameters your model actually
exposes, and adapt the `make_<Name>_hyperpar` / `<Name>_hyperpar` functions in
`R/<name>.R`.

### On strategy inheritance (issue #274)

[Issue #274](https://github.com/traitecoevo/plant/issues/274) asked for three
things. Two are **handled by the scaffolder**:

1. an explicit environment template argument (`environment =`), and
3. no redundant bindings/files/tests for reused components — a reuse-environment
   model generates **no** `<Name>_Environment` anywhere.

The third ask — having the new strategy **inherit** an existing strategy's C++
implementation (the old `FF16r : public FF16_Strategy` idea) rather than cloning
it — is **not** automated, and the scaffolder does not generate an inheriting
subclass. It clones the template's strategy code (renamed), which you then edit.
This is deliberate: `FF16r` was removed from the codebase, no strategy currently
subclasses another (all derive directly from `Strategy<E>`), and the base
methods are mostly non-`virtual`, so genuine override-inheritance would first
require making the relevant `*_Strategy` methods `virtual` in the C++ core. If
you want that pattern, treat it as a C++ design change to the base strategy
first, then hand-write the subclass — the scaffolder's job ends at wiring the
`<Strategy, Environment>` pair through the build.

## Step 3 — build, then test

```sh
make rebuild        # clean → RcppR6 → full_compile → roxygen
```

`make rebuild` is required because the yml/header changes are interface changes.
Then:

```r
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
p <- scm_base_parameters("XX24") |>
  add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
res <- run_scm(p, collect = TRUE, refine_schedule = TRUE)
# res$species is a tidy tibble (time, node, height, mortality, fecundity, …)
# Strategy parameters live under $pars, e.g. p$strategies[[1]]$pars$lma
```

A freshly scaffolded clone should reproduce its template exactly until you
change the biology — a useful sanity check. Run `devtools::test()` (or
`testthat::test_file("tests/testthat/test-strategy-<name>.R")`) and flesh out
the generated test file.

## Verification checklist

- `make rebuild` compiles and links cleanly (new `<name>_*.o` appear).
- `environment_type("XX24")` returns the expected env (`"XX24_Env"`, or the
  reused env in reuse mode).
- In reuse mode, `grep -r "<Name>_Environment"` finds **nothing** — no phantom
  environment was generated.
- `run_scm(..., collect = TRUE, refine_schedule = TRUE)` runs end-to-end.
