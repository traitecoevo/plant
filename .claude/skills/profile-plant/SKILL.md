---
name: profile-plant
description: >-
  Profile and benchmark the plant C++/R package: find hotspots, measure whether
  a change or branch helped or regressed, and compare performance across builds.
  Use when asked to profile, benchmark, time, or optimise plant; to find what is
  slow; to compare speed between two branches or commits; or to show a change
  did not regress timing. Covers building each arm into its own library, the
  interleaved A/B method (wall-clock drifts more than most changes are worth),
  instruction counts, compiled-code comparison, and native /usr/bin/sample
  hotspot analysis.
---

# Profiling and benchmarking plant

plant is a two-layer package (a templated C++ core under `inst/include/plant/`, bridged to R by RcppR6). Almost all run time is inside the C++ `.Call` boundary, so profiling is about measuring the compiled binary carefully and localising C++ hotspots. What counts as evidence of a speed change is set out in [.claude/measurements.md](../../measurements.md); this skill is the workflow.

## Rules

1. **Benchmark an installed build, one library per arm.** `R CMD INSTALL --no-docs -l <lib> .` builds with the package's own flags; `devtools::load_all()` compiles with different ones and its numbers are not comparable. Two builds of plant cannot load into one R session, so each arm is its own library and its own `Rscript`.
2. **Interleave, and compare medians.** Wall-clock on a laptop drifts by 10 to 30% between rounds of the same binary. Alternate base and head invocations for several rounds and compare each arm's median. Never compare against a number measured in another session.
3. **Read the instruction count beside the time.** `/usr/bin/time -l Rscript ...` reports instructions retired, which repeats to within about 0.3% for a fixed script where the time does not. A change that saves time must also save instructions.
4. **Know what the compiled code did.** Compare base and head objects function by function (the `compare-compiled-code` skill, with its positive control). A refactor that claims no change shows identical functions; a speed change shows exactly which functions moved.
5. **The baselines do not move.** A speed change leaves `default_baseline.rds`, `reference-gradient.tsv`, `FF16_reference/` and the scenario gateway's scorecard bit-identical unless a moved number was agreed in advance. Reordering floating-point arithmetic is a moved number.
6. **Native sampling perturbs timing.** Sample to find where time goes, never during an A/B.

## The benchmark scripts

| Script | What it times |
|---|---|
| `scripts/bench_run.R <model> [reps] [lifetime]` | One stand run on the default schedule, FF16 or TF24. A repeat is ten FF16 runs, or one TF24 run. |
| `scripts/bench_gradient.R <model> [reps] [lifetime]` | The forward run that records the states, and the reverse sweep over every metric and column, FF16 or TF24. |
| `run_plant_benchmarks()` (exported) | FF16's run and schedule refinement through `bench::mark`. |

Each prints one tab-separated line per repeat, so rounds from two builds can be pasted together and summarised in R.

## Recipe: A/B two commits

```sh
S=<scratch dir>
git worktree add --detach $S/base origin/develop
for arm in base head; do
  src=$S/base; [ $arm = head ] && src=.
  mkdir -p $S/lib-$arm
  R_MAKEVARS_USER= R CMD INSTALL --no-docs -l $S/lib-$arm $src
done
for round in 1 2 3 4 5; do
  for arm in base head; do
    R_LIBS=$S/lib-$arm Rscript scripts/bench_run.R TF24 3 | sed "s/^/$arm\t$round\t/"
  done
done > $S/times.tsv
for arm in base head; do
  R_LIBS=$S/lib-$arm /usr/bin/time -l Rscript scripts/bench_run.R TF24 3 2>&1 |
    grep 'instructions retired'
done
```

Run the gradient script the same way when the change reaches the reverse sweep or anything phylloptim computes.

## Finding hotspots

Start a long run (`bench_run.R TF24 50`), find its process id, and sample it:

```sh
/usr/bin/sample <pid> 10 -file $S/profile.txt
```

The "Sort by top of stack, same collapsed" section lists self-time symbols with sample counts; diff it between two builds to see which symbol moved. If `sample` cannot attach, it says so and the run is unperturbed. Rprof only shows that the time is inside `.Call`.

## Where the time goes

**FF16:** `SCM::run` → the stepper → `Patch::compute_rates` → `Species`/`Node::compute_rates` → `FF16_Model::compute_rates` → `net_mass_production_dt` → `assimilation`, the crown quadrature. The recurring hotspots are the quadrature integrand, the light spline lookup, competition, and libm `pow`/`exp`. What has paid off there, and is easy to undo by "tidying" (agents.md §12):

1. The quadrature integrand passed as its own closure type, not a `std::function`, so it inlines at every point.
2. Integer `eta` in the crown shape as multiplication chains rather than `pow`, and `1/height` cached in an aux slot.
3. The O(1) spline span lookup in odelia's `hermite_spline`.
4. Small hot helpers defined in headers, since the build has no link-time optimisation and an out-of-line helper is a real call per node.
5. Integer state and aux indices, never a `std::map<std::string, int>` lookup, on the hot path.

**TF24:** dominated by the leaf's hydraulic solve in phylloptim (the root-collar solve over the layered root water balance), not the crown quadrature. Improvements to the solve land in phylloptim and are timed here on a TF24 run and a TF24 sweep. To iterate quickly, shorten `lifetime`: the cost per step is unchanged, only the number of steps.

## Reporting

Put the evidence in the PR's first comment: the commands; each arm's median and the spread of its rounds; the instruction counts; the compiled-code comparison (which functions changed, and the positive control); and which baselines were checked bit-identical.
