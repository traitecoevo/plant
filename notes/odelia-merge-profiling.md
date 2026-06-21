# Profiling gate before merging `odelia` → `develop`

**Status (2026-06-22):** the `odelia` branch (plant depends on the external
`odelia` package for the ODE solver / interpolator instead of the in-tree
solver) is **~17–20% slower than `develop` on FF16 SCM**. The cause is
**not** the solver swap — it is the spline interpolator. Do not merge to
`develop` until this gap is closed (or consciously accepted).

## Root cause

`odelia::spline::Spline::operator()` does a `std::lower_bound` **binary search on
every evaluation**. develop's in-tree `tk::spline` was already optimised to an
O(1) index lookup (uniform fast path + non-uniform proportional-guess+nudge —
plant#435). The FF16 assimilation quadrature hits the environment spline at
every quadrature point × node × ODE step, so that per-lookup search dominates.

Native sample (`build_schedule`, identical 8 s window), top-of-stack self counts:

| symbol | develop | odelia (merged) | Δ |
|---|---:|---:|---:|
| spline `operator()` | `tk::spline` 1622 | `odelia::spline::Spline` 2259 | **+637 (~+39%)** |
| `ode::Step::step` (solver itself) | 21 | 19 | ~equal |

The solver step profiles equal — the regression is entirely the spline.

## Upstream fix (recoverable pointers)

- **odelia#21** — issue (commented with this evidence).
- **odelia PR #22** (`spline-uniform-index`) — implements O(1) index for uniform
  and non-uniform grids, bit-identical. **This is the fix.**
- Once #22 (and odelia#26, the `.onLoad` DLL fix) land and odelia is reinstalled,
  plant picks it up with no plant-side change. Then revert plant's workaround
  (plant#491) and re-run the A/B below.

## Baseline numbers (2026-06-22, this machine)

Same-machine, interleaved A/B, **sample off**, 7-rep one-iteration median +
20-repeat wall-clock totals. Both builds share identical FF16 strategy code.

| FF16 metric | develop | odelia (merged) | ratio |
|---|---:|---:|---:|
| one-iter `scm` (ms) | 68 | 80–81 | ~1.18× |
| one-iter `build_schedule` (ms) | 204–208 | 246–247 | ~1.20× |
| 20-rep `scm` (s) | 1.39–1.40 | 1.62–1.64 | ~1.17× |
| 20-rep `build_schedule` (s) | 4.17–4.22 | 4.93–4.97 | ~1.18× |

**Merge gate:** after odelia#22 lands + reinstall, the spline self-sample delta
and the whole-run gap should largely vanish — expect odelia(merged) within
run-to-run noise (~10%) of develop. If it does not, re-profile before merging.

## How to re-check (same-machine A/B)

Methodology follows `notes/profile-ff16-2026-06-16.md`: **sample off** for
comparing source/dependency changes (machine drift makes cross-session absolute
numbers untrustworthy — always measure both states back-to-back on the same
box). You cannot load two builds of `plant` in one R session, so build each and
run the harness in a separate `Rscript`.

```sh
# 1. develop baseline in a worktree
git worktree add -f /private/tmp/plant-develop develop
( cd /private/tmp/plant-develop && make compile )

# 2. this branch (rebuild after reinstalling the fixed odelia)
make compile

# 3. run the harness against each build, interleaved (script below)
for r in 1 2; do
  Rscript --no-init-file /tmp/bench_ab.R "$(pwd)" "odelia-r$r"
  Rscript --no-init-file /tmp/bench_ab.R /private/tmp/plant-develop "develop-r$r"
done
```

For the native breakdown (which symbol moved), use the repo's profiler with
`/usr/bin/sample` on (needs permission to inspect the R process) and diff the
"Sort by top of stack" section of the `.sample.txt`:

```sh
PLANT_PROFILE_REPEATS=30 PLANT_SAMPLE_SECONDS=8 Rscript scripts/profile-benchmarks.R FF16
```

### Harness (`/tmp/bench_ab.R`)

```r
# Same-machine A/B: FF16 scm + build_schedule (refine) cases.
# Mirrors scripts/profile-benchmarks.R. One-iteration median (ms, 7 reps) and
# 20-repeat wall-clock total (s). Sample off.
args  <- commandArgs(trailingOnly = TRUE)
path  <- if (length(args))      args[[1]] else "."
label <- if (length(args) > 1)  args[[2]] else path
suppressMessages(pkgload::load_all(path, compile = FALSE, quiet = TRUE))

mk_scm <- function() expand_parameters(trait_matrix(0.0825, "lma"),
                                       scm_base_parameters("FF16"))
mk_bs  <- function() { p <- scm_base_parameters("FF16")
                       p$strategies <- list(FF16_Strategy()); p$birth_rate <- 0.1; p }
run_scm_case <- function() { run_scm(mk_scm()); invisible(NULL) }
run_bs_case  <- function() { run_scm(mk_bs(), refine_schedule = TRUE); invisible(NULL) }

med_ms    <- function(f, n = 7) { ts <- numeric(n)
  for (i in seq_len(n)) ts[i] <- system.time(f())[["elapsed"]]; median(ts) * 1000 }
total20_s <- function(f) system.time(for (i in 1:20) f())[["elapsed"]]

invisible(run_scm_case()); invisible(run_bs_case())  # warm up
cat(sprintf("RESULT|%s|scm_one_ms=%.1f|bs_one_ms=%.1f|scm_20_s=%.3f|bs_20_s=%.3f\n",
            label, med_ms(run_scm_case), med_ms(run_bs_case),
            total20_s(run_scm_case), total20_s(run_bs_case)))
```

## Note on the dependency build

This branch links plant against odelia's compiled DLL for the XAD `Tape` runtime
(`src/Makevars` + macOS `install_name_tool` fixup, and `R/zzz.R`). That is a
stopgap pending odelia#26 (plant#491); it does not affect the timings above
(the link/load is one-time, not in the hot path).
