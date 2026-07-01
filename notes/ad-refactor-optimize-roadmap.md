# AD gradient machinery: refactor, consolidate, optimize

## Context

Reverse-mode automatic differentiation of emergent trait gradients is **built and
validated end-to-end** for FF16, TF24, and TF24f on branch `spike-ff16-scm-emergent`
(101 commits, ~21,700 net lines over `develop`). The prototype phase is done; this plan
is the **refactor + optimize phase**. It does **not** add modelling capability — it
makes the validated machinery faster, smaller, more correct, and CI-covered.

Three problems motivate it:

1. **The R-side `Rcpp::as<>` env round-trip is the correctness ceiling.** The gradient
   is two-pass: pass 1 runs a resident `run_scm(save_RK45_cache=TRUE)` and harvests a
   frozen schedule + per-RK-stage environment; pass 2 is an AD replay that never re-runs
   the solver. Pass-1 harvest happens **in R** (`ff16_harvest`, `tf24f_harvest`, …),
   pulling `patch$environment_history` across the boundary and rebuilding each env via
   `Rcpp::as<>`. That reconstruction is **not faithful** for the crown-sampled light
   above cohort heights — it is the root cause of the TF24f census-recon fidelity floor
   at patch lifetime > 4 and of a validation trap (a native-correct AD `g'` read as "18%
   wrong" only through a round-tripped env). It also rebuilds the whole RcppR6 patch on
   every access (O(stand size), ~1600× slower per `R/emergent_gradient.R:70-73`).

2. **Duplication.** Three near-identical replay engines (`src/ff16_emergent.cpp` 1851,
   `src/tf24_emergent.cpp` 479, `src/tf24f_emergent.cpp` 2266) plus three R harvest
   seams. The Cash-Karp stepper and state containers are *already* shared; what is
   triplicated is the **orchestration glue** (the `Frozen` struct, `build_frozen`/`as<>`
   marshalling, per-strategy deriv/axpy/seed, `census_reduce`, coupled-canopy
   reconstruction).

3. **No benchmarks or CI coverage for the AD path.** There is no timing harness for the
   gradient calls, and 7 FF16 AD test files `skip()` with "AD tape symbols unavailable
   in this load_all session" — so the FF16 AD surface is not exercised in CI. ~31
   intermediate `ad_*`/`coupled_*`/`recon_*` validation scripts accumulated during the
   prototype.

**Intended outcome:** the AD path drives native environments (round-trip gone, the
long-horizon fidelity floor dissolved), one strategy-parameterized engine instead of
three, a benchmark + regression net that proves every step preserves the validated
Jacobians, and the spike merged to `develop` as a clean reviewable PR stack.

### Decision on PR #541 / merge sequencing (resolved with Dan)

**Do NOT gate this work on landing #541 first; close #541 without merging.** #541
(`spike-ff16-hierarchy`, the `<T,E,S>` templating) is already the ancestor base of this
spike — nothing depends on it being in `develop`, and its commits remain in the spike
history whether or not the PR is open. The only value of keeping it open was reviewer
ergonomics, which is better served by re-cutting fresh PRs from the *final* tree.
**Action: close PR #541 (do not merge)** — it isn't serving a purpose as a live review
object. We do all the work in place on the spike (no rebase, no merge conflict), then
decompose the final tree into a stacked PR sequence for `develop` (Phase 5), where a
foundation slice (the templating) can be cut fresh as PR #1. *(Closing the PR is a
mutating GitHub action — done on plan approval, not during planning.)*

---

## Phase 1 — Safety net: benchmarks + correctness-regression fixture (do first)

Establish the baseline *before* touching anything, so every later step is provably
non-regressing in both timing and value. Templates the existing
`scripts/bench_ab.R` / `scripts/bench_tf24.R` conventions.

- **`scripts/bench_gradient.R`** — same-machine A/B harness (mirrors `bench_ab.R`:
  `path`+`label` args, `pkgload::load_all(compile=FALSE)`, warm-up, one `RESULT|…` line
  per build, interleave this-branch vs a worktree). The novel requirement: **time the
  three sub-costs separately** by calling internals directly —
  - `run_ms` = resident `run_scm(save_RK45_cache=TRUE)` (context),
  - `harvest_ms` = the R-side `*_harvest()` call alone (**the refactor's target**),
  - `impl_ms` = the C++ `*_impl` replay+sweep alone,
  - headline `harvest_frac = harvest_ms/(harvest_ms+impl_ms)` and a `val=<digest>`
    bit-identity check (à la `bench_tf24.R:31`).
  Canonical cases lifted from test fixtures (fixed schedules): `ff16_frozen` (28 traits,
  4 metrics), `ff16_resident_coupled`, `ff16_resident_ms`, `ff16_offspring`,
  `tf24f_census`, `tf24f_resident`, `grow_individual`. `BENCH_REPS`/`BENCH_CASES` env
  vars; optional `BENCH_SCALING=1` to factor `impl_ms ≈ replay(traits) + M·sweep`.

- **Correctness-regression fixture (AD-vs-AD baseline).** Snapshot the current validated
  Jacobians (values + `d(metric)/d(trait)` matrices) to
  `tests/testthat/fixtures/gradient-baseline.rds` for every distinct engine
  (frozen / coupled-single / coupled-ms / state-jacobian / offspring / grow, across
  FF16/TF24/TF24f). A new `tests/testthat/test-gradient-regression.R` reloads it and
  asserts two-tier tolerance: **bit-identical** (`1e-12` rel) for pure-relocation steps,
  **noise-floor** (`~5e-6` rel) for the coupled/ms paths that legitimately reorder FP
  sums. This is complementary to the existing AD-vs-FD tests (those pin AD to physics at
  ~1%; this pins AD to its own validated self at machine precision). Use fixed node
  schedules so snapshot/check are reproducible.

- **Extend `scripts/profile-benchmarks.R`** with gradient cases (`grad_ff16_frozen`,
  `grad_ff16_resident`, `grad_tf24f_census`) reusing its Rprof + native `/usr/bin/sample`
  machinery. Per the profile-plant skill: sample **ON** for hotspot localisation, **OFF**
  for A/B ratios; only same-session ratios are trustworthy.

- **Timing history:** recorded alongside the `bench_gradient.R` RESULT lines, columns
  `date | step+sha | case | run_ms | harvest_ms | impl_ms | public_ms | harvest_frac |
  cum_speedup | incr_speedup | sample | fixture(PASS/FAIL + worst rel_dev) | notes`.
  No speedup row is recorded without its correctness verdict attached.

**Validates:** the fixture itself is the validation infrastructure; confirm it passes on
the current HEAD and that `harvest_ms + impl_ms ≈ public_ms`.

---

## Phase 2 — Native-env harvest in C++ (kills the `Rcpp::as<>` round-trip)

> **UPDATE 2026-06-30 (measured, supersedes the fidelity rationale below).** The premise
> that the `Rcpp::as<>` env round-trip causes the TF24f census floor is **NOT supported by
> measurement**. A lossy-vs-native A/B of the census recon gives bit-identical LAI at every
> horizon (H=4 rel 6.3e-7; H=5 3.98e-3; H=8 1.19e-2, no abort) — `native == lossy` exactly.
> The TF24_Environment round-trip is faithful for the census path (knots fully determine the
> spline; `compute_competition(0)` is read at z=0, inside the domain). **The lifetime>4 floor
> is the g' backward-FD near-cancellation + frozen-schedule replay, not the env.** Phases 2a
> (FF16) and 2b (TF24f census) are DONE and committed as **bit-identical refactors** that
> remove a real serialization round-trip + the documented env-validation hazard and unlock the
> `harvest_frac` perf headroom — **not** a fidelity fix. Lifting the census floor is a separate
> item (exact-AD g' / refined-schedule replay). No false lifetime gate added. See memory
> `ad-env-roundtrip-faithful-for-census`.

The highest-value change: it has a **correctness payoff** (removes the lifetime-floor and
the validation trap) and shrinks the per-engine surface Phase 3 then unifies. The
round-trip is structurally avoidable — confirmed: the replay reads the env only as a
*double* value + analytic derivative (`src/ff16_emergent.cpp:125-128`), lifting just those
scalars into the AD type; the AD type never touches the env. A faithful native env
**pointer** is sufficient and exact.

**Approach (engine stays OUTSIDE the SCM object).** Do **not** add a `collect_gradient`
run mode to `SCM::run()` — the harvest is
already captured during the normal `save_RK45_cache=TRUE` run. Instead:

- Add `inst/include/plant/gradient/resident_harvest.h` with a `ResidentHarvest` struct
  and a free `harvest_resident(const Patch<T,E,S>&, int species)` builder that gathers,
  **by const-ref/pointer into the Patch's own storage**, the pieces R rebuilds today:
  `step_history`, `&environment_history[n][s]` (native pointers, not copies),
  `stand_*_stage_history`, `stand_newnode_*`, and per-species cohort metadata. Move the
  two remaining R computations — trapezoid node weights (`emergent_gradient.R:94-97`) and
  the per-RK-stage `pr_survival` matrix (`:99-103`) — into this C++ builder.
- Change each engine's `Frozen.eh` from owning `std::vector<…Environment>` to holding
  `const Environment*` borrowed from the Patch. **Delete the `Rcpp::as<>` loop**
  (`ff16_emergent.cpp:450`, and the TF24/TF24f equivalents). This line is the round-trip.
- The C++ entry takes the **live `Patch&`** from the RcppR6 external pointer (env pointers
  are only valid while the SCM C++ object lives; RcppR6 keeps it alive in R). The R API
  (`stand_gradient()` etc.) becomes a thin forwarder passing the live SCM xptr down.

**Sub-steps, each bit-checkable against Phase 1's fixture:**
1. **2a** — `ResidentHarvest` + C++ gather + the two moved arithmetic computations, still
   building `Frozen` from copied envs. **Bit-identical** to current `ff16_harvest`
   outputs (assert `tw`, `ppsurv`, sampled env values equal). Isolates the gather move.
2. **2b** — switch `Frozen.eh` to borrowed native pointers; delete `as<>`. Bit-identical
   at lifetime ≤ 4; at **lifetime 5–8** the recon LAI must now match
   `compute_competition(0)` to ~1e-6 and TF24f must no longer abort at the spline edge —
   add new long-horizon gates to `test-tf24f-census-gradient.R` that were previously
   impossible. **This is the headline correctness win.**
3. **2c** — repeat for TF24 / TF24f env entries.

**Files:** new `inst/include/plant/gradient/resident_harvest.h`; `inst/include/plant/patch.h`
(public accessors already exist, `:153`, `:176-218`); `src/ff16_emergent.cpp`,
`src/tf24_emergent.cpp`, `src/tf24f_emergent.cpp` (`Frozen`, `build_frozen`, `attach_*`);
`R/emergent_gradient.R`, `R/tf24_emergent_gradient.R`, `R/tf24f_emergent_gradient.R`
(harvests shrink to forwarders). `make rebuild` if an export signature changes.

**Validates:** all `test-{ff16,tf24,tf24f}-*gradient*.R` stay green; the new lifetime-5/8
TF24f gates pass; the Phase-1 fixture confirms ≤4-horizon Jacobians unchanged; record the
`harvest_frac` drop in the timing-history file.

---

## Phase 3 — Unify the three replay engines (maximal dedup, per Dan)

> **OUTCOME 2026-06-30 (safe dedup done; finding on the limit of further unification).**
> Executed the safe, high-value, bit-identical shares and measured the result. What was
> genuinely duplicated and is now shared in `inst/include/plant/gradient/`:
> - `scm_harvest.h`: `recover_birth_rate` (was inlined 5×), `birth_steps` (2×),
>   `census_trapezium` (the descending-height census reduction, FF16 + TF24f).
> - `coupled_canopy.h`: `canopy_comp_at` (the Yokozawa light trapezium, was hand-copied
>   as ff16 `coupled_comp_at` + tf24f `tf24f_comp_at`; 4 call sites).
> - The Cash-Karp stepper (`ff16_cashkarp_replay`) + state containers were **already**
>   shared before this work.
>
> All bit-identical (fixture max_rel 0.0; one 9e-16 FP-reassociation on frozen census).
> **Key finding:** the engines shrank only ~107 lines for ~132 lines of shared header —
> i.e. the duplicated *orchestration* was small relative to the **irreducible
> strategy-specific physiology** (FF16's light-response hyperbola + deep-crown GK vs
> TF24f's tracked-collar leaf solve + curvature harvest; the coupled `deep_net_coupled`
> / deriv kernels). So the "4900 → 2000" estimate was optimistic: after sharing the
> stepper/harvest/census/canopy, the bulk of each engine is **distinct biology, not
> redundant copy**. A generic `replay_cohort_final<Strategy,S>` trait-class would add
> machinery to abstract over genuinely-different deriv kernels for modest line savings,
> and the multi-species coupled unification is high-risk (the stiffness/conditioning
> guards) for low dedup value (the MS engines are largely distinct). **Recommendation:
> stop the engine unification here** — the redundant duplication is removed; what remains
> is either irreducible or net-negative to force. The remaining `build_interp`
> single-vs-MS intra-file overlap is a small (~25-line) optional same-file cleanup.
> Original (pre-execution) plan kept below for context.

Collapse the triplicated orchestration into one strategy-parameterized engine. The
stepper (`ff16_cashkarp_replay`) and state containers are already shared across all three
(`tf24_emergent.cpp` and `tf24f_emergent.cpp` already `#include ff16_production_kernel.h`).
The engines differ in `replay_cohort_final` at exactly three points: the env→net call, the
state-vector width (5/6/7), and the rate fill.

**Minimal trait-class interface** each strategy supplies (mostly *binding* things already
templated — `FF16ProdPars<S>`/`TF24ProdPars<S>` and `*_compute_rates_from_net` already
exist with matching shapes):
```
struct StrategyTraits<Strat,S> {
  using ProdPars  = …;            // existing
  using LifeState = …;            // strategy-supplied width (TF24f adds collar + log_density)
  static S        net_at(ProdPars, Env*, LifeState, /*harvest@n,stage*/);
  static Rates    rates_from_net(ProdPars, height, area_leaf, net);
  static LifeState seed(ProdPars, Env* birth_env, h0, …);   // establishment / IFT collar birth
};
```
TF24f's collar-as-6th-state + curvature harvest fits as a **hook**, not a fork: the collar
is an extra `LifeState` component with rate `k_acclim·dprofit_dpsi`; the curvature harvest
is a *harvest producer* orthogonal to the replay template. **Do not** merge TF24's
envelope-zeroed collar with TF24f's taped collar — they are different math; keep them as
two `seed`/`net_at` implementations behind the hook.

**Ordered dedup (each gated by Phase-1 fixture); Dan chose to push through all ranks:**
1. **Shared `census_reduce` + `(w,f)` metric dispatch** — pure arithmetic, currently
   hand-copied identically (`ff16:587`, `tf24f:226`). **Bit-identical.** (~300 lines)
2. **Generic `replay_cohort_final<Strat,S>`** via the trait class (frozen/invasion path).
   FF16/TF24 offspring are bit-exact references. **Bit-identical.** (~400-600 lines)
3. **Generic coupled-canopy reconstruction** — the Yokozawa trapezium
   `Q=(1-(z/h)^eta)^2` is identical across strategies (`coupled_comp_at`/`build_interp`
   vs `tf24f_comp_at`/`build_canopy`); only the per-cohort rate differs. Validated to the
   ~1e-4 FD-noise floor + sign-flip invariants (not bit-identical — the coupled tape is
   numerically delicate). (~500-800 lines)
4. **Multi-species engine** (`assemble_metrics_coupled_ms`/`FrozenMS`). **Highest risk** —
   the MS stiffness/conditioning guards (joint-env-drift > 1e-3 gate, finiteness guard,
   tape scoping) must be carried through **verbatim**; a refactor that drops them silently
   re-opens divergence. Validate against the MS-R0/MS-cross-species-R1 gates on fast fixed
   schedules.

**Files:** `inst/include/plant/models/ff16_production_kernel.h` (shared `census_reduce`,
generic `replay_cohort_final`, trait base), new
`inst/include/plant/gradient/replay_engine.h` (strategy-agnostic orchestration); the three
`src/*_emergent.cpp` shrink to: include the generic engine + supply the strategy trait
class + the `[[Rcpp::export]]` entry points.

**Validates:** ranks 1-2 bit-for-bit; rank 3 at ~1e-4 + sign-flips; rank 4 against the MS
gates. Run the **whole** test suite each step (generic tests loop over
`helper-plant.R` strategy lists; remember TF24f is the non-5-state variant excluded
there). Record line-count reduction + timing per step.

---

## Phase 4 — Cleanup: fold validation scripts into tests, cover the AD surface in CI

Per Dan: **fold the still-valuable validation scripts into proper tests**, then delete the
scripts (rather than deleting outright or leaving them as loose scripts).

- **Migrate** the meaningful checks from `scripts/coupled_r0_check.R`,
  `coupled_r1_check.R`, `recon_static_check.R`, `coupled_vs_full_scm_fd.R`,
  `coupled_multispp_r0.R`, `coupled_multispp_r1.R`, `coupled_gap_decompose.R`,
  `coupled_birthenv_channel.R`, etc. into `testthat` tests (or fold into the existing
  `test-{ff16,tf24f}-*gradient*.R`). Delete the migrated scripts. Keep a small curated set
  of `ad_*` demos as runnable examples; the `overstorey-staging` guide `.qmd` remains the
  user-facing narrative.
- **Make the AD tests actually run in CI (enabler — otherwise folded tests just skip).**
  The 7 FF16 AD test files skip with "AD tape symbols unavailable in this load_all
  session"; TF24f tests already run live. Diagnose the asymmetry and make the FF16 AD
  symbols available under the test load path — tie to the **odelia DLL load-ordering**
  fix (plant needs `importFrom(odelia,…)` so odelia's `.onLoad` runs before plant's
  `useDynLib`) and the build-optimized-then-`load_all(compile=FALSE)` path. The aim:
  folded + existing AD tests execute in CI rather than skip. (Aligns with roadmap Phase C:
  "make the headline result CI-testable in plain R".)

**Validates:** the migrated tests pass under `make test` / `devtools::test()` *without*
skipping; the AD surface is genuinely covered. Confirm no orphaned `source(scripts/…)`
references remain.

---

## Phase 5 — Decompose for merge to `develop`

PR #541 is **closed (not merged)** at the start of this work — its commits stay as spike
ancestors. Re-cut the final spike tree into a stacked, reviewable PR sequence (this is
*merge* decomposition, independent of the in-place writing order above):

1. **Foundation** (`<T,E,S>` templating + exact AD growth-rate gradient) — additive,
   bit-identical to current package; cut fresh from the final tree, the clean foundation
   reviewers vet first (replaces the role the old #541 would have played).
2. **Native-env harvest** (Phase 2) — the correctness fix; land even if later phases slip.
3. **Engine unification** (Phase 3 ranks 1-2 bit-identical; rank 3 noise-gated; rank 4 MS).
4. **Cleanup + CI coverage** (Phase 4).

Each PR carries its slice of the timing-history table and its fixture verdict.

---

## Verification (end-to-end)

- **Build discipline:** C++-only change → `make compile`; new export / RcppR6 field →
  `make rebuild`; then `pkgload::load_all(".", compile=FALSE)` (optimized -O2, **never**
  `devtools::load_all()` which is -O0). Run the **whole** plant suite before each PR
  (`make test`), not just `test-strategy-*` — generic tests loop over `helper-plant.R`.
- **Correctness gate (every step):** `Rscript scripts/gradient_fixture.R check` (or the
  `test-gradient-regression.R` testthat run) — bit-identical for relocation/rank-1-2
  steps, noise-floor for coupled/ms. Plus the existing AD-vs-FD tests stay green.
- **The Phase-2 headline proof:** new lifetime-5 and lifetime-8 TF24f census gates that
  were impossible under the round-trip now pass to ~1e-6 vs `compute_competition(0)`.
- **Timing (same-session A/B):** build this branch + a worktree of the pre-step ref,
  interleave `Rscript scripts/bench_gradient.R <path> <label>` runs, compare `RESULT`
  lines; expect `harvest_frac` → ~0 and `public_ms` to drop by the old harvest fraction on
  the frozen/offspring/grow cases. Native `/usr/bin/sample` (OFF for ratios) confirms time
  left the R harvest. Record cumulative + incremental speedup alongside the bench results.

## Critical files

- `R/emergent_gradient.R` (`ff16_harvest` :61, `ff16_harvest_ms` :129, the O(stand) comment
  :70-73, `stand_gradient` :229) — the R-side round-trip + arithmetic being moved to C++.
- `R/tf24_emergent_gradient.R`, `R/tf24f_emergent_gradient.R` — TF24/TF24f harvests + FD gates.
- `src/ff16_emergent.cpp` (reference engine: `Frozen` :66-117, `build_frozen` :442-462,
  `replay_cohort_final` :225-256, coupled :643-921), `src/tf24_emergent.cpp`,
  `src/tf24f_emergent.cpp` (7-state collar replay :71, curvature harvest, coupled canopy).
- `inst/include/plant/patch.h` (native `environment_history` :153, `stand_*_stage_history`
  :176-218 — the harvest source).
- `inst/include/plant/models/ff16_production_kernel.h` (shared stepper/state — where the
  generic engine + trait class live); `tf24_production_kernel.h`.
- New: `inst/include/plant/gradient/resident_harvest.h`, `…/gradient/replay_engine.h`.
- New: `scripts/bench_gradient.R`, `tests/testthat/fixtures/gradient-baseline.rds`,
  `tests/testthat/test-gradient-regression.R`; extend `scripts/profile-benchmarks.R`.

---

## Jobs for the NEXT session (handoff 2026-06-30, session 2)

Session 2 landed (branch `spike-ff16-scm-emergent`, 5 commits on top of `2156a42f`;
suite FAIL=0 ERROR=0 SKIP=9 **PASS=2587**, fixture all PASS):
- `8153c924` — **the old Job 1 NaN, fixed.** Root cause was NOT a `log_density` runaway —
  it was trait **`a_l2`**: the last cohort's birth step lands on the final step
  (`birth==N`), which the coupled `rn<N` loop never established, leaving a zero-height
  cohort whose `area_leaf=(h/a_l1)^(1/a_l2)` differentiates to `0*log(0)=NaN` (and biased
  the value). Fixed by establishing `birth>=N` cohorts at `h0` (FF16 + TF24f). All 28 FF16
  traits finite, single + MS + **>2 species**; re-snapshotted 3 fixtures.
- `279cae5c`, `cde17580` — **the old Job 2 (dead-code), done.** Removed 7 superseded
  `_impl`/probe entries; reframed `bench_gradient.R` to public/native only; purged 21
  `ad_*.R` sourceCpp prototypes; guide made self-contained.
- `bef2503e` — **FF16 multi-species coupled gradient natived** (`ff16_coupled_gradient_ms_native`
  + `build_frozen_ms_scm`; public path off the R `ff16_harvest_ms`; `_impl` kept as the
  test FD-reference; bit-identical, `_impl`/`_native` share `ff16_coupled_gradient_ms_core`).
- `b6830bbe` — **notes/ pruned to just this roadmap**; all `see notes/X` pointers trimmed.

> **SESSION 3 DONE (2026-06-30).** Tasks A and B below are both complete; 3 commits on
> `49a163a1` (suite FAIL=0 SKIP=9 **PASS=2591**, fixture all PASS).
> - **A done (`e82841f2`):** TF24f MS coupled gradient natived. Extracted
>   `tf24f_coupled_gradient_ms_core` (C++ structures) from the monolith; `_impl` (R-list,
>   the FD-reference) is a thin wrapper; new `tf24f_coupled_gradient_ms_native` reads the
>   live Patch; wired `tf24f_resident_census_gradient_ms_ad` to it. native == impl
>   bit-for-bit (diff 0.0). Cosmetic/perf, no result change.
> - **B done (`47e24bad`):** single-species coupled long-horizon NaN root-caused + gated.
>   Both handoff suspects RULED OUT by probe (collar curvature k_acclim*d2p_dpsi2 ~ -9.9,
>   stably decaying; g' light channel — freezing it doesn't help, breaks H4). Real cause:
>   the coupled **log_density<->canopy sensitivity feedback** goes unstable — deeply-shaded
>   TF24f cohorts carry a LARGE leaf `dprofit_dL` (7->24) + high density, amplifying the
>   canopy-reshaping feedback; step-cap bisection localises ignition to a single step
>   (jac -2.1 -> +4e9, compounding to 1e+237/NaN with horizon). The SCM tames it via ADAPTIVE
>   stepping; the frozen-step replay cannot (double env_err ~1e-6@H4 -> ~2e-2@H5+). FF16
>   robust (closed-form net, bounded light response). Fix: the single-species path had NO
>   gate; added the double R0 (env_err) gate + post-sweep finiteness/magnitude guard -> clear
>   error for stiff horizons. H<=4 unchanged. True horizon extension = adaptive sub-stepping
>   in the replay (future hardening).
>   - **Plausibly related to #550** (`[TF24 hydraulics] SCM cohort-density blow-up under
>     extreme seasonal drought`): SAME size-density equation `d(log_density)/dt = -dg/dh -
>     mortality` destabilised by a large/near-singular `dg/dh` originating in the TF24 leaf.
>     But a DISTINCT manifestation, confirmed by measurement: #550 is a genuine `run_scm`
>     **caustic** (the density VALUE overflows to +Inf, log_density ~1e9-1e10, under drought;
>     Dan: "not a stepper artifact", no gradient exists there). Ours has `run_scm` HEALTHY in
>     the failing regime (steady env, max|log_density| ~5, nowhere near overflow) and the TRUE
>     resident gradient EXISTS (full-SCM FD finite at H=6: d(LAI)/d(lma)=16.9) — only the
>     frozen-step AD *replay* diverges. So ours is a replay artifact of the same equation, not
>     the caustic; #550's root fixes (NSC #517 / hydraulic-optimum continuity) would also help
>     ours, and #550's regime is a hard limit for gradients too (no derivative at a caustic).
>     Note also TF24f uses the smooth tracked collar, avoiding the #550 opt_psi_stem
>     discontinuity — so ours is the milder cousin.
> - **C (birth_rate deriv) — FF16 resident slice DONE (Dan chose "build FF16 resident slice
>   now").** `birth_rate_gradient(scm, metrics, species)` (exported) returns the
>   resident-coupled `d(census metric)/d(birth_rate)` for FF16 (LAI/biomass/size_moment).
>   Mechanism: `assemble_metrics_coupled` templated on the birth_rate type; new
>   `ff16_birth_rate_gradient_core/_native` register birth_rate as the SOLE tape input on the
>   same coupled whole-stand replay the trait gradient uses (one reverse sweep/metric). The
>   FROZEN reading is the trivial identity `metric/birth_rate` (verified op/br=0.86; offspring
>   stays frozen even under resident -> excluded); only the resident axis needs a tape, and it
>   is genuinely non-trivial -- the canopy feedback dominates and FLIPS the sign of biomass
>   (AD biomass -6.25e-4 vs frozen identity +0.456). Validated AD == FD over the SAME coupled
>   recon (perturb birth_rate arg of `ff16_coupled_metrics_impl`) to ~0.7%. New fixture case
>   `ff16_birth_rate` (noise tier) + test in test-ff16-stand-gradient.R.
> - **C cross-species block DONE** (`8840d874`): `assemble_metrics_coupled_ms` gains an
>   optional active per-species birth_rate vector; `ff16_birth_rate_gradient_ms_core/_native`
>   register only the target species' birth_rate; `birth_rate_gradient(scm, metrics, species)`
>   dispatches single/multi + gates the MS path on env_err. The cross term
>   `d(total-stand metric)/d(birth_rate_s)` is the genuinely new object (frozen = diagonal
>   identity). AD == FD over the same coupled MS recon ~0.3%; MS test added.
> - **C demographic-equilibrium axis DONE** (`51ca9616`): exact-AD
>   `d(net_reproduction_ratio)/d(birth_rate)` = dR0/db, via Dan's MUTANT framing (mutant
>   traits = resident; change in mutant fitness as resident density moves = the density
>   feedback). `assemble_metrics_coupled` carries a survival-weighted offspring accumulator
>   (CensusState -> FullState, decoupled => census bit-identical); new metric
>   "net_reproduction_ratio" with UNIT per-seed weights tw/birth_rate0 (birth_rate0 = fixed
>   harvest br, new Frozen field, so birth_rate enters R0 only via the canopy). Single-species
>   only (cross-species census-only). R0 recon == SCM exactly; dR0/db < 0; AD==FD(same recon,
>   weight held) 0.19%; ~25% to full-SCM FD = frozen-grid response. Suite PASS=2610; test-ad
>   PASS=362. The plant-side access point for the regnans R0=1 Newton solve is complete.
>   **Remaining (future):** TF24f resident birth_rate (gated like its trait gradient); the
>   generic stand_gradient metric-loop wiring; cross-species net_reproduction_ratio if needed.

**Two tasks chosen for the next session (Dan, 2026-06-30):**

### A. Native the TF24f multi-species coupled gradient (mirror `bef2503e`)

The lone public path still on the R harvest. `tf24f_resident_census_gradient_ms_ad`
(`R/tf24f_emergent_gradient.R`) calls `tf24f_harvest_ms` (R, rebuilds the RcppR6 patch
O(stand) + `Rcpp::as<>` env) then `tf24f_coupled_gradient_ms_impl` twice (gate_only=TRUE
then FALSE). **Cosmetic/perf only — it works and is tested; changes no results.** Pattern
to copy from FF16 (`src/ff16_emergent.cpp`): (1) extract `tf24f_coupled_gradient_ms_core`
from `tf24f_coupled_gradient_ms_impl` (`src/tf24f_emergent.cpp:1956`) — but note that impl
is a ~400-line MONOLITH that builds `EH`/`NNH`/`birth` inline from R lists (`:2005-2040`),
unlike FF16 where `build_frozen_ms` was already a separate fn; so the refactor is to move
that R-list→C++ build into the `_impl` wrapper and have the core take the C++ structures;
(2) add a native builder reading `patch.environment_history` / `step_history` /
`stand_newnode_*_stage_history_all` / `plant::gradient::birth_steps(patch,s)` (the live
`Patch`); (3) add `tf24f_coupled_gradient_ms_native(SEXP scm_, pp_list, ...)` that bundles
the R0 gate (env_err) + the AD sweep like the FF16 native does; (4) wire
`tf24f_resident_census_gradient_ms_ad` to it. **`compileAttributes()` + `make compile`**
(new export). Gate: native == impl bit-for-bit on the MS fixture stand; fixture +
whole suite stay green.

### B. TF24f long-horizon coupled-gradient NaN (a REAL correctness gap)

`tf24f_resident_census_gradient_ad` (single-species) goes **all-NaN** on patch lifetime
**H≥6** (H=4 fixture is fine — it has `birth==N` too, so this is NOT the boundary-cohort
bug A above, which is already fixed). Diagnosis gathered this session:
- The **frozen census recon is clean at H8** (`log_density` ∈ [−2.8, 5.6], all heights >0),
  so the NaN is purely in the **COUPLED AD tape**, and it is **derivative-only** (the
  `values` are finite; only the adjoint is NaN — a single op with finite value but
  ∞/NaN local derivative poisons the whole shared tape).
- It **scales with horizon** (H4 ok, H6/H8 NaN). **FF16 coupled is robust to lifetime 140**,
  so it is TF24f-leaf-solve-specific, not generic.
- Prime suspects: the `g' = (height_dt - g_back)/GEPS` backward-FD near-cancellation
  accumulated over more stages (the handoff's original `log_density`/g' hypothesis, but
  for the COUPLED tape), or a `pow`/`anchor`/interpolator-derivative singularity in
  `build_canopy` (`src/tf24f_emergent.cpp`) at a knot for the taller H8 stand.
- Method: instrument the AD core (`tf24f_coupled_gradient_core`) to find the first op whose
  derivative goes non-finite (e.g. guard + report `xad::value`/`xad::derivative` per stage,
  or bisect by zeroing channels). Compare H4 (works) vs H6 (fails) at the same code path.
  Consider a `log_density`-rate cap mirroring `check_initial_density_rates`, or reading the
  canopy derivative more carefully where a knot coincides with a cohort crown edge.

### C. (still pending) Derivative w.r.t. birth_rate — see §3 below (NEW capability).

---

## Jobs from session 1 (handoff 2026-06-30) — §1 DONE, §2 DONE; §3 pending

Recorded at the end of the refactor/optimize + API-consolidation session. Phases 1–4
are done and committed; the suite is green and the AD-vs-AD fixture is bit-identical.

### 1. Multi-species coupled resident gradient — make it reliable

> **DONE 2026-06-30 (FF16 fully; TF24f boundary-cohort fixed, one TF24f follow-up left).**
> Root cause of the default-28-trait NaN was **not** a `log_density` runaway — it was a
> single trait, **`a_l2`**, on LAI/biomass: a cohort whose introduction time lands on the
> final step (`birth_step == N`, which always happens for the last node at
> `max_patch_lifetime`) was **never established** by the coupled re-evolution loop
> (`for rn < N`), so its census height stayed at the zero-initialised value. The frozen
> engine's per-cohort `replay_cohort_full` establishes such a cohort at `h0` then runs zero
> Cash-Karp steps, so it never had a zero-height cohort — that asymmetry was the bug. With
> a zero height, `area_leaf = (h/a_l1)^(1/a_l2)` is differentiated at `h=0`, giving
> `0*log(0) = NaN` in the `a_l2` column (FF16 localised it there; TF24f's shared census
> tape was poisoned to all-NaN). The zero-height cohort also **biased the census VALUE**
> (the coupled total only matched the frozen-summed reference to ~5e-3 before; ~1e-7 after).
> **Fix:** after the `rn < N` loop, establish any still-unborn cohort (`birth >= N`) at `h0`
> in the frozen final birth env, never stepped — mirrors the frozen engine exactly.
> Applied to FF16 (`assemble_metrics_coupled`, `assemble_metrics_coupled_ms`) and TF24f
> (all four coupled loops: `tf24f_coupled_metrics_impl`, `tf24f_coupled_gradient_core`, and
> the two MS loops in `tf24f_coupled_gradient_ms_impl`). Validated: all 28 FF16 traits
> finite single + multi-species; **> 2 species** works (3-species, every target, R0 value
> == frozen-summed to ~1e-8); `a_l2` AD == FD-of-recon to ~1e-3. New regression tests in
> `test-ff16-stand-gradient.R` (all-traits-finite + >2-species) and
> `test-tf24f-census-gradient.R` (boundary-cohort). Re-snapshotted `ff16_resident`,
> `ff16_resident_ms`, `tf24f_resident` (the old snapshots had the corrupted-value bug).
> Full suite FAIL=0 ERROR=0 SKIP=9 PASS=2587; fixture all PASS. FF16 robust to long stands
> (lifetime 140). **Refined-schedule stiffness gate unchanged** (still fixed-schedule-only).
>
> **REMAINING TF24f follow-up (separate, deeper, NOT the boundary cohort):** TF24f's
> *single-species* coupled (resident) census **gradient** still goes all-NaN on **longer
> horizons** (H≥6; H=4 fixture works — it has `birth==N` too and is fine, proving this is
> not the boundary cohort). The frozen census recon is clean at H8 (`log_density` well-
> behaved, heights all > 0), so the NaN lives in the **coupled AD tape** specifically — a
> derivative-only NaN (the VALUE is finite) that scales with horizon. Prime suspects: the
> `g' = (height_dt - g_back)/GEPS` backward-FD near-cancellation accumulated over more
> stages, or a `pow`/`anchor`/interpolator-derivative singularity in `build_canopy` at a
> knot for the taller H8 stand. This is the TF24f equivalent of the numerically-delicate
> coupled tape the Phase-3 note flagged; needs its own diagnosis (instrument the AD core to
> find the poisoning op). The FF16 coupled tape does **not** have this (robust to lifetime
> 140), so it is TF24f-leaf-solve-specific.

`stand_gradient(feedback = "resident")` on a multi-species stand (the cross-species
total `d(total-stand metric)/d(θ_s)`) is **not fully functional**:

- It works for a **small trait set** on a 2-species **fixed** node schedule (validated
  with `traits = c("lma", "a_p1")` in `test-ff16-stand-gradient.R` and the guide's
  multi-species chunk), but the **default 28-trait** sweep drives some cross-species
  feedbacks to **NaN** — so it is currently pinned to a hand-picked trait subset.
- **> 2 species is untested.**
- On an adaptively-**refined** schedule it can go stiff (a cohort's `log_density` runs
  away); the public path gates this with a clear error (see the "instability" section
  of the guide), so it is fixed-schedule-only.

Goal: get it running reliably **for all traits** and **for > 2 species**. Likely work:
diagnose which traits NaN and why (the `log_density` runaway / `g'` near-cancellation in
the joint re-evolution is the prime suspect), consider the `log_density`-rate cap
mirroring the SCM's `check_initial_density_rates`, and add `> 2`-species test coverage.

### 2. Dead-code audit — flag anything not used by the core interface

> **AUDIT DONE 2026-06-30 (caller-grep over every `[[Rcpp::export]]` in the 3 emergent
> TUs, against `R/` non-generated + `tests/` + `scripts/` + C++ internal callers). The
> handoff's "orphaned `_impl`, no callers" framing is only partly right — 3 of the 6
> prime candidates are LIVE in the bench harness, and the `_impl` entries are the
> perturbable-by-`pp` FD-reference surface the validation tests need (the `_native` twins
> read a live Patch and can't be re-evaluated at a perturbed trait on the frozen harvest).
> Verdict:**
> - **REMOVED (commit follows this note): 7 entries** — the 4 truly-dead
>   (`ff16_census_reconstruct_impl`, `tf24_stand_gradient_impl`, `tf24f_census_recon_impl`,
>   `ff16_reverse_tape_probe`) + the 3 that were used **only** by `bench_gradient.R`'s
>   `impl_ms` timing path (`ff16_coupled_gradient_impl`, `tf24f_census_gradient_ad_impl`,
>   `tf24f_coupled_gradient_impl`). The bench was **reframed to call the public (native)
>   entries only** (`public_ms` + the `val` regression digest; the migration-era
>   `harvest_ms`/`impl_ms` split is gone — native fuses the harvest into the C++ call, so
>   there is no separable R harvest left to time). `compileAttributes()` + `make compile`
>   regenerated the bindings; fixture all PASS, suite FAIL=0 ERROR=0 SKIP=9 PASS=2587.
> - **KEPT — the multi-species resident path's SOLE implementation (NOT legacy — never
>   natived):** `ff16_coupled_gradient_ms_impl`, `ff16_coupled_metrics_ms_impl`,
>   `tf24f_coupled_gradient_ms_impl` (the public `stand_gradient(feedback="resident")` MS
>   path calls these directly). Removing them removes the MS gradient; eliminating their
>   R harvest needs a native MS entry (a real follow-up, the migration single-species got).
> - **KEPT — FD-reference surface in tests** (perturb `pp` on a frozen harvest, the
>   "AD == FD over the SAME reconstruction" checks): `tf24f_offspring_gradient_impl`
>   (`test-tf24f-census-gradient.R:153`), `ff16_state_jacobian_impl`,
>   `tf24_state_jacobian_impl`, `ff16_coupled_metrics_impl`, `tf24f_coupled_metrics_impl`,
>   `ff16_coupled_metrics_ms_impl`, `tf24_offspring_production_gradient_impl`.
> - **`_core` (C++-internal, exported but R=0 by design):** `tf24f_coupled_gradient_core`,
>   `tf24f_offspring_gradient_core` — called by their `_impl` + `_native`; keep (could
>   drop the `[[Rcpp::export]]` to de-register the unused R binding, a minor follow-up).
>
> So the "collapse `_impl`/`_native`/`_core` triples" is NOT a clean win: the `_impl`
> (raw-`pp`) and `_native` (live-Patch) entries serve genuinely different callers (FD
> validation + benchmarking vs the fast public path). Recommended action: remove only the
> 4 truly-dead entries; keep the rest as the deliberate test/bench surface.

The native-harvest migration left a layer of `_impl` (R-list, `Rcpp::as<>`-env)
back-compat entries that the public R API no longer calls (the API now routes to the
`_native` entries). Audit all the new AD code and flag for likely removal. Concrete
starting points (each "no R/test caller" per a `grep` of `R/` + `tests/` — **verify**
before deleting, some are still used by FD-reference tests):

- Orphaned exported `*_impl` C++ entries (prime candidates):
  `tf24_stand_gradient_impl` (TF24 offspring now routes to `*_native`),
  `ff16_coupled_gradient_impl`, `ff16_census_reconstruct_impl`,
  `tf24f_census_gradient_ad_impl`, `tf24f_census_recon_impl`,
  `tf24f_coupled_gradient_impl`. (Keep the `_native` + the `_core` they share; remove
  the dead `_impl` wrappers + their RcppR6 exports.)
- The `_impl` / `_native` / `_core` triplication generally: once an `_impl` has no
  caller, collapse to `_core` + `_native`.
- `scripts/ad_*.R` (≈21 sourceCpp prototypes) + the kept demo scripts: superseded by the
  compiled paths + the guide. **Decision (Dan):** scripts/examples may be **kept for the
  commit onto `develop`** (so they are in history) and **deleted in a follow-up** if we
  judge the history value worth it — i.e. don't delete pre-merge; revisit post-merge.
- The 7 FF16 sourceCpp AD test files that `skip()` under load_all/CI ("AD tape symbols
  unavailable") — decide: convert to exercise the compiled path, or retire (the compiled
  path is already covered by `test-ff16-stand-gradient` + the regression fixture).

Method: a `grep`-for-callers pass over every `[[Rcpp::export]]` and new R function;
anything reachable only from removed/▾dead paths is a removal candidate. The Phase-1
fixture + full suite gate every deletion (must stay FAIL=0, bit-identical).

### 3. Derivative w.r.t. birth_rate of all species (new capability)

A new gradient *axis* (not a trait): `d(emergent metric)/d(birth_rate_s)` for every
species' birth-rate driver. **Use case (Dan):** evolving the community toward
**demographic equilibrium** — birth_rate enters the offspring weighting (`tw ∝
birth_rate`) and the establishment initial condition, so its derivative gives the
Newton/gradient step that drives the stand to self-replacement (e.g.
`net_reproduction_ratio → 1`, or a target `offspring_production`). Additive to the
existing tape: register the per-species `birth_rate` as an AD input alongside (or
instead of) the traits in the reverse sweep — the harvest already carries it as a
scalar. Likely a `birth_rate = TRUE`-style option on `stand_gradient` returning the
metrics × species birth-rate sensitivities (and the cross-species block for the
multi-species coupled path).
