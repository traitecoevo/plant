---
name: plant-update-interface
description: >-
  Migrate code that uses the plant R package up to plant's latest interface.
  Use when a downstream product (analysis scripts, a dependent package, a
  notebook) breaks or needs updating after plant changed — e.g. "update X to the
  latest plant", "plant's interface changed, fix my code", "migrate to the new
  run_scm/Control API", "bring this repo up to date with plant". Reads plant's
  NEWS.md "Breaking changes" as the spec, finds affected call sites, and
  proposes per-site fixes. Does NOT change plant itself or invent API mappings —
  if NEWS lacks a needed mapping, it stops and asks.
---

# Updating products to plant's latest interface

plant's R-facing interface changes as the package develops (the SCM
consolidation, the `Control()` defaults, the odelia split, …). Downstream
code — analysis repos, dependent packages, notebooks — breaks when symbols are
renamed, removed, or change meaning. This skill migrates that downstream code.

**The contract this skill works under:**

- **`NEWS.md` is the spec, not this file.** plant's `NEWS.md` "Breaking changes"
  section holds the canonical `old -> new` mappings, written to be machine-read.
  This skill is the *applier*: it reads NEWS each run, so it keeps working as
  plant adds changes without this file having to be updated. **Never migrate
  from a mapping memorised here or from training data** — always read the live
  NEWS. If a breakage has no mapping in NEWS, stop and ask the maintainer rather
  than guessing.
- **Propose, don't blind-apply.** Show the diff for each call site. Mechanical
  renames can be applied directly; semantic changes (an argument's meaning
  flipped, a default changed) must be flagged for human review.
- **This runs in the *product* repo, not in plant.** It edits the downstream
  code. It only *reads* plant (its NEWS, optionally its installed source).

## Step 1 — locate the authoritative NEWS.md

Resolve plant's NEWS in this order; stop at the first that works:

1. **Working inside the plant repo itself?** (`./DESCRIPTION` has
   `Package: plant`) → use `./NEWS.md`.
2. **Latest from GitHub** (the usual case — migrating a separate product to the
   newest plant):
   ```sh
   curl -fsSL https://raw.githubusercontent.com/traitecoevo/plant/develop/NEWS.md
   ```
   Pin to a tag/branch the product is targeting if one is known
   (`.../plant/<ref>/NEWS.md`).
3. **Installed package** (offline / matching what's actually loaded):
   ```r
   system.file("NEWS.md", package = "plant")   # or: utils::news(package = "plant")
   ```

Read the whole **"Plant (development version)"** block — at minimum its
"Breaking changes" subsection, plus "New features" (a product may want to adopt
one, e.g. `export_patch_state()` or `control$shading_model`).

## Step 2 — determine the product's baseline

You need to know which changes are *new to this product*, so you don't re-apply
ones already done.

- **Marker file** (preferred): look for `.plant-interface-version` at the product
  root. The skill writes it at the end of a successful run, recording the date
  and the list of breaking-change keys already applied. If present, only entries
  **not** in it are in scope.
- **No marker?** Fall back to the installed `packageVersion("plant")` and/or ask
  the maintainer "which plant version did this last work against?". The plant
  dev version doesn't bump per-change, so version numbers are coarse — when in
  doubt, treat the *entire* current "Breaking changes" section as candidate work
  and let Step 4's grep tell you what actually applies (a no-match symbol is
  simply already-clean).

## Step 3 — build the migration map from NEWS

From the in-scope "Breaking changes" entries, extract each `old -> new` mapping.
Classify each:

- **Mechanical rename** — a 1:1 symbol or simple-arg substitution
  (e.g. `fast_control()` -> `control()`). Safe to apply across all matches.
- **Semantic change** — same symbol, changed meaning/default, or a multi-step
  replacement (e.g. raw `Control()` now means *fast* not *accurate*; the
  `run_scm()` argument-order change; `run_scm_error()` -> set a field then read
  another). These need per-site judgement — never sed blindly.

If an entry's mapping is ambiguous or missing, **stop and ask** — do not invent
one. That gap is a NEWS bug; note it so it can be fixed at the source.

## Step 4 — find the call sites

Grep the product for every affected symbol. Search R sources, scripts, vignettes
/ Quarto / Rmd, tests, and roxygen examples:

```sh
grep -rn --include='*.R' --include='*.r' --include='*.Rmd' --include='*.qmd' \
  -e 'run_scm_collect' -e 'run_scm_error' -e 'build_schedule' -e 'split_times' \
  -e 'fast_control' -e 'scm_base_control' -e '\brun_scm\b' -e '\bControl\b' \
  -e 'OdeRunner' .
```

**Known affected symbols as of 2026-06** (operational grep list — *reconcile
against the live NEWS Breaking changes section, which is canonical; add any
symbol it lists that isn't here, and don't trust this list if NEWS is newer*):

| Symbol to grep | Kind | See NEWS entry |
|---|---|---|
| `run_scm_collect`, `run_scm_error`, `build_schedule`, `split_times` | removed → `run_scm()` flags | SCM consolidation (#408/#459/#462) |
| `run_scm(` (positional args) | arg-order change | SCM consolidation |
| `fast_control`, `scm_base_control` | removed → `control()` / `control_accurate()` | Control defaults (#463) |
| `Control(` / `control(` used expecting *accurate* tolerances | semantic: now fast | Control defaults (#463) |
| `OdeRunner`, plant ODE/interpolator C++ headers | moved to odelia | odelia split (#456/#464) |

Also grep `NAMESPACE` / `DESCRIPTION` `Imports`/`importFrom` if the product is a
package, and any `library(plant)` / `plant::` qualified calls.

## Step 5 — apply the migrations

Per call site, smallest faithful edit, in the product's own style:

- **Mechanical renames:** apply directly. Show a summary diff.
- **Semantic changes:** present the before/after for the maintainer. Examples:
  - A bare `Control()` / `scm_base_control()` that fed a *high-accuracy* run must
    become `control_accurate()`, not `control()` — you cannot tell which was
    intended from the call alone, so ask or flag.
  - `run_scm(p, env, ctrl, use_ode_times = TRUE)` written positionally may now
    bind `use_ode_times` to the wrong parameter — rewrite with named arguments.
  - `run_scm_error()` becomes a two-step (`scm$collect_errors <- TRUE`; read
    `scm$combined_node_errors`) — restructure, don't substitute.
- Leave a brief `# plant N.N migration: …` comment only where the change is
  non-obvious; don't litter.

## Step 6 — rebuild, test, record

Run the **product's own** verification (not plant's):

- Package → `devtools::load_all()` then `devtools::test()` (or `R CMD check`).
- Scripts / notebooks → execute the affected ones end-to-end.
- Re-grep the removed symbols to confirm none remain (except in NEWS/changelog
  prose).

Then update the baseline marker so the next run is incremental:

```r
# .plant-interface-version (product root) — written after a clean migration
writeLines(c(
  "migrated_on: 2026-06-23",
  "plant_ref: develop@<sha or tag>",
  "applied_breaking_changes:",
  "  - scm-consolidation",
  "  - control-defaults",
  "  - odelia-split"
), ".plant-interface-version")
```

Report what changed, what was applied mechanically vs. flagged for review, and
anything that couldn't be migrated because NEWS lacked a mapping.

## Verification checklist

- Migration mappings came from the **live NEWS**, not memory — and every applied
  change traces to a NEWS "Breaking changes" entry.
- No removed symbol remains in the product's executable code.
- Semantic changes (accurate-vs-fast `Control`, `run_scm` arg order) were
  surfaced for review, not silently rewritten.
- The product's own tests / scripts pass after migration.
- `.plant-interface-version` records the new baseline.
- Any missing/ambiguous NEWS mapping was reported (so plant's NEWS can be fixed).
