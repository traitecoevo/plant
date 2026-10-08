---
name: profile-plant
description: >-
  Profile and benchmark the plant C++/R package — find hotspots, measure whether
  a change or branch helped or regressed, and compare performance across builds.
  Use when asked to profile, benchmark, time, or optimize plant; to find what is
  slow; to compare speed between two branches/commits; or to verify a change did
  not regress timing. Covers the build-then-load workflow, the sample-off
  same-machine A/B methodology (machine drift is large — only same-session ratios
  are trustworthy), native /usr/bin/sample hotspot analysis (Rprof alone is too
  coarse for C++), and the timing-history reporting convention.
---

# Profiling & benchmarking plant

plant is a two-layer package (templated C++ core under `inst/include/plant/`,
bridged to R via RcppR6). Almost all run time is inside the C++ `.Call`
boundary, so profiling is really about **measuring the compiled binary
carefully** and **localising C++ hotspots**.

## Golden rules (read first — these are the mistakes that waste hours)

1. **Build, then load the built DLL — never benchmark a `load_all()` rebuild.**
   ```sh
   make compile        # builds src/plant.so with the package's -O2 flags
   ```
   Then load with `pkgload::load_all(compile = FALSE, quiet = TRUE)` (the
   scripts below do this). `devtools::load_all()` with recompilation can use
   different flags and gives non-comparable numbers.

2. **Comparing a change → run with native sampling OFF.** `/usr/bin/sample`
   runs concurrently with `Rprof()` and inflates the Rprof totals, which has
   repeatedly produced false "regressions". Use `PLANT_PROFILE_REPEATS=20`
   *without* `PLANT_SAMPLE_SECONDS` for A/B. Turn sampling ON only when you need
   the native call stack to find *where* the time goes.

3. **Machine drift is large; only same-session, same-machine ratios are
   trustworthy.** Absolute ms/s wander 10–30% between sessions. Always measure
   both states **back-to-back, interleaved**, recompiling between them. Never
   compare today's numbers to numbers written in a note last week — re-measure
   the baseline now.

4. **Verify the change is bit-identical (or know why it isn't).** Run targeted
   strategy/individual tests after each change. If you reorder floating-point
   arithmetic (e.g. reciprocal-multiply), say so and check whether any test
   tolerance had to move. Bit-identical refactors are the safe kind.

5. **You cannot load two builds of `plant` in one R session.** For cross-branch
   work, build each in its own checkout/worktree and run one `Rscript` per
   build.

## The two canonical FF16 cases

Both live in `scripts/profile-benchmarks.R` and `scripts/bench_ab.R`:

- **`scm`** — `run_scm(p)` for an FF16 resident. The core SCM solve.
- **`build_schedule`** — `run_scm(p, refine_schedule = TRUE)`. Adds the C++
  cohort-refinement loop; assimilation-quadrature heavy, so it is the more
  sensitive case for solver/spline/strategy changes.

## Recipes

### A. Profile one build (find current hotspots)

```sh
make compile
# Rprof totals + bench::mark, sample off (the comparable mode):
PLANT_PROFILE_REPEATS=20 Rscript scripts/profile-benchmarks.R FF16
# Add native call stacks (may need permission to inspect the R process):
PLANT_PROFILE_REPEATS=30 PLANT_SAMPLE_SECONDS=8 Rscript scripts/profile-benchmarks.R FF16
```
Outputs go to `tmp/profile-benchmarks/<timestamp>/`:
`benchmark-summary.csv`, `*-by-self.csv` / `*-by-total.csv` (Rprof),
`*.sample.txt` (native).

### B. A/B a source change (same branch, before/after)

```sh
git stash            # or commit the change
make compile && Rscript scripts/bench_ab.R . before
git stash pop
make compile && Rscript scripts/bench_ab.R . after
```
Run each 2–3× and compare the medians; ignore differences inside ~10% noise.

### C. A/B across branches (e.g. this branch vs develop)

```sh
git worktree add -f /private/tmp/plant-other develop
( cd /private/tmp/plant-other && make compile )
make compile                                   # this branch
for r in 1 2; do
  Rscript --no-init-file scripts/bench_ab.R "$(pwd)"               "this-r$r"
  Rscript --no-init-file scripts/bench_ab.R /private/tmp/plant-other "other-r$r"
done
```
Interleaving controls for machine drift. `git worktree remove /private/tmp/plant-other` when done.

### D. TF24 single-species timing

`scripts/bench_tf24.R` times one MPL=10 TF24 SCM run and prints a trajectory
metric (so you can confirm bit-identity alongside timing). `BENCH_REPS` sets the
rep count.

## Reading the results

- **Rprof is too coarse for C++.** `*-by-self.csv` will show ~97–99% in
  `.Call` / `SCM___FF16__FF16_Env__run` — that only confirms time is in the core.
- **Native `/usr/bin/sample` is where the real signal is.** In `*.sample.txt`,
  the **`Sort by top of stack, same collapsed`** section lists leaf (self-time)
  symbols with sample counts — diff this section between two builds to see which
  symbol moved. Same sampling-window length on both → counts are comparable.
- If `sample` prints `sample cannot examine process …; try running with sudo`,
  it failed to attach; the Rprof run is then *un*-perturbed (which is fine for
  A/B), you just don't get native stacks.
- **FF16 hot path:** `SCM::run` → `ode::Solver::step` → `derivs` →
  `Patch/Species/Node::compute_rates` → `Individual::growth_rate_given_height` →
  `FF16_Strategy::compute_rates` → `net_mass_production_dt` → `assimilation`
  (the crown quadrature). The recurring hotspots are the assimilation quadrature
  callable, the environment **spline lookup** (`operator()`), competition, and
  libm `pow`/`exp`.

## Where FF16 time goes (levers that have worked)

A 2026-06 optimisation pass took FF16 ~2–3.5× faster end-to-end. The targets,
hottest first, with the techniques that paid off — start here before
re-deriving:

1. **Assimilation-quadrature callable.** `QK::integrate` is a template; pass the
   integrand as a bare lambda (its closure type), not a `std::function`, so it
   inlines at every quadrature point. Biggest single FF16 win.
2. **Crown-shape `pow`.** `Q(z,H)` raises `(z/H)` to integer `eta`; specialise
   common `eta` (1,2,4,8,10,12) as multiply chains (`CanopyShape`, issue #465)
   instead of libm `pow` (issue #361), and hoist/cache the `1/height` reciprocal
   (in aux) instead of recomputing it per quadrature point.
3. **Environment spline lookup.** Hit per quadrature point × node × step. Use the
   O(1) index — uniform fast path plus non-uniform proportional-guess+nudge —
   instead of `std::lower_bound` (issue #435). This is exactly what odelia#21 /
   odelia PR #22 port into the odelia package that the `odelia` branch depends on.
4. **Inline hot helpers across TU boundaries.** The no-LTO build forces a
   cross-TU call for every out-of-line helper; moving small FF16 hot helpers
   (`area_leaf`, `compute_competition*`) into the header lets the compiler fold
   them into the per-node loops. (LTO itself was tried — issue #470 — for no
   measurable gain, and declined.)
5. **No `std::map<string,int>::at` in hot state/rate paths.** Resolve aux/ode
   indices once; use direct integer indices.
6. ~~**Gradient calc.** Reuse a thread-local scratch `Individual`.~~ **Withdrawn.**
   It measured within 1.5% of the plain copy, and `growth_rate_gradient` is now
   scalar-templated, so a cached scratch would hold active scalars and carry tape
   slots between recordings. Copy per call; see `plant/agents.md`.

Most of these were bit-identical (or a deliberate, documented reciprocal-multiply
reorder).

**TF24 has a different hot path.** It is dominated by the plant hydraulic solver
(`Leaf::find_root_collar_psi`, a Brent inner loop over the multi-layer root water
balance), not the assimilation quadrature. Profile it with `scripts/bench_tf24.R`;
to iterate quickly drop `max_patch_lifetime` (per-step cost is unchanged, only the
number of steps). `set_physiology` is *not* on the hot path.

## Reporting convention

For each change report: the exact command; the one-iteration
`scm`/`build_schedule` medians and the 20-repeat totals; **cumulative** speedup
(vs the original baseline) and **incremental** speedup (vs the previous state);
whether native sampling was on or off; and the bit-identity check (which tests
passed). Keep a running timing-history table — in the PR or a fresh `notes/`
file — so each change's gain is attributable, not just the aggregate.

The full prior optimisation log (per-change timings, the sampler-effect audit,
per-target detail) was `notes/profile-ff16-2026-06-16.md`, removed so this skill
is the single source of truth. Recover it from git history:
`git show a2e5e2ff:notes/profile-ff16-2026-06-16.md` (or
github.com/traitecoevo/plant/blob/a2e5e2ff13c87b6561ea07ac412d0ada71d112b0/notes/profile-ff16-2026-06-16.md).

## Need a full flamegraph?

`/usr/bin/sample` (above) is the primary native path and is enough for almost
all hotspot work. If you need a full call-graph image, the usual heavier tools
apply — Xcode Instruments (macOS), Linux `perf` / `uProf`, or `gperftools` /
`jointprof` + `pprof` with Rcpp. None is wired up in the repo; set them up
against current upstream docs when needed.

## Files

| path | what |
|---|---|
| `scripts/profile-benchmarks.R` | canonical driver: bench::mark + Rprof (+ optional native sample) for FF16/K93/TF24 |
| `scripts/bench_ab.R` | same-machine A/B harness (one-iter median + 20-rep total) for two builds |
| `scripts/bench_tf24.R` | TF24 single-species timing + trajectory metric |

## Environment variables (driver)

| var | default | meaning |
|---|---|---|
| `PLANT_PROFILE_REPEATS` | 1 | repeats inside each profiled case (use 20+ for stable totals) |
| `PLANT_SAMPLE_SECONDS` | 0 | seconds of `/usr/bin/sample`; 0 = off (use off for A/B) |
| `PLANT_RPROF_INTERVAL` | 0.005 | Rprof sampling interval (s) |
| `PLANT_BENCHMARK_ITERATIONS` | 1 | `bench::mark` iterations in `run_plant_benchmarks` |
| `BENCH_REPS` | 5 | timed reps in `bench_tf24.R` |
