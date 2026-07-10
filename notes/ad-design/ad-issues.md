# AD infrastructure: work items, dependencies, and build order

Companion to [`ad-infrastructure-design.md`](./ad-infrastructure-design.md).
Each item is scoped to one layer (odelia / plant / UX) or is a de-risking
prototype. Items are deliberately small; where a design is still open the item
says what must be decided and by which prototype.

**Legend.** Class: **CP** critical path · **PROTO** de-risking prototype (gates a
CP item) · **NTH** nice-to-have. Repo: `odelia` / `plant` / `meta` (docs, fixtures).

---

## Prototypes (do first — each can move a decision)

### PROTO-1 — Scalar-cost measurement · repo: plant · gates PLANT-1, decision 2
Implement the FF16 frozen cohort as a uniform-`value_type` odelia System, gradient
one metric, compare `tape.getMemory()` and wall-clock to the spike's mixed engine
on the canonical cases. **Output:** a number that chooses uniform scalar vs. an
encapsulated "frozen input" fallback. Blocks committing to the full-scalar Patch.

### PROTO-2 — TF24 census cross-sensitivity · repo: plant · gates TF24 census, decision 7
One TF24 cohort; gradient a density-dependent census metric with the leaf-optimizer
IFT delivered as an on-tape `SuppliedDerivative`; check against finite differences.
**Output:** confirmation (or refutation) that the full-scalar `run_mutant` replay
captures the density→optimum cross-term. **Highest risk in the whole design** — run
early; a refutation rescopes TF24 census out of v1.

### PROTO-3 — `SuppliedDerivative` / `CheckpointCallback` ergonomics · repo: odelia → plant · feeds ODELIA-3
Minimal `CheckpointCallback` that injects a known d(out)/d(in) into a reverse tape;
reproduce one leaf-optimizer sensitivity from the spike bit-for-bit. **Output:** the
`SuppliedDerivative` API shape. Prerequisite for the clean version of PROTO-2.

---

## odelia layer (generic, plant-agnostic)

> **Status: this whole layer is implemented** (ODELIA-1..6, RIF-1..3, PROTO-3), on a
> stack of PRs on the odelia fork. The items below are kept for the dependency map and
> build order; the settled names are `DifferentiationTargets` (was `Independents`),
> `SuppliedDerivative` (was `AnalyticEdge`), `record_stage`/`record_ode_step`/`replay_step`
> + `has_recorded_field()` (was `cache_*`/`has_cache`), and the `Replayable` concept.
> Follow-ups live in the odelia tracker: **#22** (interpolator unification), **#23**
> (history rows), **#27** (functional = pure reduction; driver owns the replay), **#28**
> (pare the demonstrator, retire `live|frozen`), **#25/#26** (comments, tests).

### ODELIA-1 — Generalize `compute_gradient` to an arbitrary functional · CP
Replace the hard-coded `sum_of_squares` with a caller-supplied functional
`std::vector<S> f(const System& replayed)`. Keep `sum_of_squares` as the prebuilt
`least_squares` functional instance. **Depends on:** —. **Blocks:** ODELIA-2, PLANT-4.

### ODELIA-2 — `compute_jacobian` (reverse, row-sweep) · CP
Record once, one adjoint sweep per output row (adjoint `XAD::computeJacobian`
pattern), reusing the Solver-owned tape. **Depends on:** ODELIA-1. **Blocks:**
PLANT-4, UX-1.

### ODELIA-3 — `DifferentiationTargets` + `SuppliedDerivative` · CP
The "differentiate w.r.t. what" bundle (`DifferentiationTargets`) and the
`SuppliedDerivative` injection built on `CheckpointCallback`. **Depends on:** PROTO-3. **Blocks:**
PLANT-2 (traits), PLANT-6 (leaf edge).

### ODELIA-4 — odelia-native test harness for the AD surface · CP-support
The spike's Jacobian test harness lives in plant; the machinery moving into odelia
needs coverage **independent of plant**, on odelia's own example systems (Lorenz,
leaf_thermal). Scope (implementation may be deferred to another dev, but the tests
are release-gating):
- **`compute_gradient` / `compute_jacobian`** — value + gradient against finite
  differences on a closed-form system; Jacobian shape/row-sweep on a multi-output
  functional; tape reuse across repeated calls. **Depends on:** ODELIA-1, ODELIA-2.
- **`SuppliedDerivative` (ODELIA-4b)** — an off-tape root-find whose IFT edge is injected
  via `CheckpointCallback`, checked against the analytic derivative and against
  differentiating the solve directly. **Depends on:** ODELIA-3.
- **Functional shape (ODELIA-4b)** — a custom functional (not `sum_of_squares`) plus
  the `least_squares` functional instance, to prove the seam is functional-
  agnostic. **Depends on:** ODELIA-1.
- **Record/replay-fixed (ODELIA-4a)** — a system with an adaptive interpolator and a
  quadrature: assert the fixed-node replay reproduces the adaptive run's value and
  that gradients match finite differences (design §7.5–7.6). **Depends on:** ODELIA-6.

These are the odelia-side equivalents of the spike's plant Jacobian fixture; plant's
AD-vs-AD oracle (UX-2) still gates the plant ports separately.

### ODELIA-5 — Document the AD API contract in `ARCHITECTURE.md` · CP-support
Extend the existing Tape-linking contract to cover the AD *API* (functional shape,
`DifferentiationTargets`, supplied derivatives) so plant depends on a versioned odelia surface. **Depends
on:** ODELIA-1..3.

### ODELIA-6 — Record-adaptive / replay-fixed numerics (the one replay primitive) · CP
Make odelia's adaptive numerics support "record node placement on the double pass,
replay on fixed nodes with the active scalar" uniformly (design §7.5–7.6). **Detailed
design: [`ad-record-replay.md`](./ad-record-replay.md)** — positions-vs-values, one
`Replayable` concept (runtime frozen/live mode, not a type split), no `Recording` noun,
and the RIF-3 anchor settled on the `Solver` member. Two mechanisms, both partly
present:
- **Record via opt-in System hooks.** The stepper calls `record_stage` /
  `record_ode_step` / `replay_step` when a System provides them (the `Replayable`
  concept, C++20 `requires` + `if constexpr`; zero-cost no-op otherwise). The payload
  is knot positions (per step) and frozen field values (per stage), not
  `stand_stage_history`.
- **Replay-fixed components.** `advance_fixed` (stepper) and `basic_interpolator<S>`
  (frozen knots, active values) exist; the **gap** is a scalar-templated fixed-rule
  `QK<S>` that consumes a recorded QAG subdivision — or a prototype showing a single
  fixed rule suffices. `QK<S>` stays single-layer (no interpolator-style wrapper).

**Depends on:** —. **Blocks:** PLANT-4a, PLANT-5a. **Tested by:** ODELIA-4a.

---

## plant layer (surgical changes to existing types)

### PLANT-1 — Uniform `value_type = S` on Patch/Species/Node/Individual · CP
Collapse the mixed active/frozen `<T,E,S>` into one uniform scalar. **Depends on:**
PROTO-1 (decides uniform vs. fallback). **Blocks:** PLANT-3, PLANT-4.

### PLANT-2 — Single scalar-templated Strategy parameter store · CP
Remove the dual double-`pars` / lifted-active-struct representation so a named
trait can be registered active. **Depends on:** —. **Blocks:** PLANT-3.

### PLANT-3 — Patch `set_params` / `set_initial_state` (odelia AD contract) · CP
Implement the System contract on Patch; map trait names + birth-rate to registered
active inputs. **Depends on:** PLANT-1, PLANT-2, ODELIA-3. **Blocks:** PLANT-4.

### PLANT-4 — Differentiate `run_mutant` (invasion gradient) · CP
The FF16 invasion gradient as `compute_jacobian` through the existing `run_mutant`
with `S=active`: canopy read **frozen** from `environment_history` (derivative
through it is zero), `is_mutant_run` suppresses self-competition. Retires
`ff16_emergent.cpp`. **Depends on:** ODELIA-2, PLANT-3. **Blocks:** PLANT-5, PLANT-7.

### PLANT-4a — Resident/total gradient: re-run the canopy on recorded knots · CP
The resident gradient on the same frozen L0/L1 schedule with the canopy **re-run
live**: re-run `compute_environment` on the *recorded* light-spline knots (§7.5)
with the active cohorts, so a trait re-shades the stand. **Do NOT** read the frozen
env (that silently yields the invasion gradient, missing self-shading) and **do NOT**
build `stand_*_stage_history` — record only the knot positions; cohort values come
from the replay. **Depends on:** PLANT-4, PLANT-5, PLANT-5a, ODELIA-6. **Blocks:**
PLANT-7 (census resident).

### PLANT-5 — Scalar-template `Species::compute_competition` + census; reuse · CP
Template the reductions on `S`; delete `gradient/{coupled_canopy.h, scm_harvest.h}`.
Enables the resident/total gradient (active canopy). **Depends on:** PLANT-1.
**Blocks:** PLANT-7.

### PLANT-5a — L2 replay: moving-node `QK` + frozen light-spline knots · CP
Wire the census height-integral to the scalar-templated Gauss–Kronrod `QK`
(`qk.h`, already written for #472) so the active plant-height bound propagates
through the *moving* quadrature nodes; reconstruct the resident light on odelia's
frozen-knot differentiable spline (design §7, L2). This is the level a frozen-node
replay misses. **Depends on:** PLANT-5. **Blocks:** census gradients (UX-3 metrics).

### PLANT-6 — Leaf-optimizer `SuppliedDerivative` for TF24/TF24f · CP
Route the forward-mode leaf sensitivity through the odelia edge. **Depends on:**
ODELIA-3, PROTO-2. **Blocks:** PLANT-7 (TF24/TF24f).

### PLANT-7 — Port TF24 / TF24f; delete parallel engines · CP
TF24/TF24f as the same differentiated `run`/`run_mutant`; remove
`tf24_emergent.cpp`, `tf24f_emergent.cpp`, the R harvest, and the five local tapes.
**Depends on:** PLANT-4, PLANT-5, PLANT-6.

### PLANT-8 — Frozen light on odelia's differentiable spline · CP-support
Build the cached resident light on `basic_spline<S>` so schedule-freeze is a
primitive. **Depends on:** PLANT-5. Folds into PLANT-4/5.

### PLANT-9 — Un-skip the 7 FF16 AD tests · NTH-then-CP
Once AD runs on odelia's tape/load path, convert the skipped tests to exercise the
compiled path. **Depends on:** PLANT-4.

### PLANT-10 — Birth-rate gradient + `d R0/d birth_rate` (equilibrium) · CP (FF16)
`birth_rate_gradient(scm, metrics, species)` on the coupled resident replay (the
frozen part is the identity `metric/birth_rate`; the resident axis needs a tape and
can flip signs), plus `d(net_reproduction_ratio)/d(birth_rate)` = dR0/db via the
mutant framing — the plant-side derivative for the R0 = 1 equilibrium Newton solve
(design §6.2). FF16 single- then multi-species in v1; TF24f gated like its trait
gradient. **Depends on:** PLANT-4a. **Blocks:** downstream equilibrium/invasion tools.

### PLANT-11 — Boundary / zero-height cohort fix + test · CP (correctness)
Establish `birth ≥ N` cohorts at the seed height `h0` so
`area_leaf = (h/a_l1)^(1/a_l2)` does not differentiate to `0·log(0) = NaN` (design
§11). Carry the spike's fix; add an explicit AD test at the final-step boundary.
**Depends on:** PLANT-4. Cheap but load-bearing — the NaN also biased the value.

---

## UX / API / workflow

### UX-1 — `stand_gradient()` stable surface · CP
Keep the public signature; forward to `compute_jacobian` + `EmergentFunctional`.
**Depends on:** ODELIA-2, PLANT-4.

### UX-2 — Regression oracle fixture · CP (do alongside migration)
Snapshot the spike's validated Jacobians to `tests/testthat/fixtures/
gradient-baseline.rds`; two-tier tolerance (bit-identity / noise floor). **Depends
on:** —. **Gates:** every PLANT-* merge.

### UX-3 — Metric set + kernels · resolved
Shipped set: **LAI, biomass, basal area** (census, reusing the scalar-templated
`compute_competition`) and **offspring_production**. Each is one `psi` kernel behind
the `EmergentFunctional` shape; mean-height and other additions are tutorial-level
extensions (design R-interface §6.4), not release scope.

### UX-4 — Gradient benchmark harness · NTH
`scripts/bench_gradient.R` timing the sub-costs; also feeds PROTO-1. **Depends
on:** —.

---

## R boundary (see [`ad-r-interface.md`](./ad-r-interface.md))

Governing invariant: only `double` crosses the R boundary; active types are
C++-internal and ephemeral. These items are the R-facing side of the driver work.

### RIF-1 — odelia `Solver_gradient`/`Solver_jacobian`/`value_and_gradient` on the double handle · CP
Retire the `bool active` flag and the R-visible `ActiveSystemType` XPtr; the driver
owns the active system internally. Include a combined `value_and_gradient(p)` that
returns both from one recording (optimizer loops call `fn`/`gr` separately; sharing
the tape halves the work — see user story §6.3). **Depends on:** ODELIA-1, ODELIA-2.
Removes the type-confusion hazard in `Solver_fit_impl`.

### RIF-2 — odelia `rebind` lift (double → active) · CP
A System contract so the driver constructs the active system generically instead of
per-example hand-construction. **Depends on:** —. **Blocks:** RIF-1.

### RIF-3 — odelia tape/active-solver cache on the double Solver · CP-support
Tape reuse within a Jacobian and across optimizer calls, without R seeing it. The
cached active solver is typed via the System's `rebind` (no `void*`/`static_cast`).
**Depends on:** RIF-1.

### RIF-4 — odelia policy: no `wrap`/`as` for active types · NTH
Keep `xad::value`/`xad::derivative` the only extraction; prevent accidental
derivative loss.

### RIF-5 — plant single `stand_gradient_cpp` entry (RcppR6 handle, C++ dispatch) · CP
One hand-written export taking the RcppR6 SCM (pointer unwrap only), dispatching
strategy/feedback in C++. **Depends on:** ODELIA-2, PLANT-4.

### RIF-6 — plant native-pointer harvest; delete `Rcpp::as<*_Environment>` · CP
The active replay reads the live Patch by pointer; no serialisation round-trip.
**Depends on:** PLANT-4, PLANT-5. (Same work as the native harvest in those items,
viewed from the boundary.)

### RIF-7 — plant thin `stand_gradient()` R wrapper · CP-support
Remove R-side branching across native/impl/resident variants. **Depends on:**
RIF-5.

---

## Build order (critical path)

```
PROTO-3 ─► ODELIA-3 ─┐
ODELIA-1 ─► ODELIA-2 ─┼─► PLANT-3 ─► PLANT-4 ─► PLANT-4a ─► PLANT-7
PROTO-1 ─► PLANT-1 ───┘        ▲          │         ▲
PLANT-2 ──────────────────────┘          ├─► PLANT-5 ─► PLANT-5a ┘
ODELIA-6 ─► (PLANT-4a, PLANT-5a) ─────────┤   PROTO-2 ─► PLANT-6 ─► PLANT-7
UX-1 ◄─ PLANT-4 ; ODELIA-4/4a/4b ◄─ ODELIA-1..3,6 (odelia-native tests)
UX-2 (oracle) ── gates every PLANT-* merge ; UX-3 resolved (input)
```

**Sequencing notes.**
- The three prototypes are independent and come first; PROTO-2 is the one that can
  rescope the release, so front-load it.
- odelia's foundation (ODELIA-1..3), the record/replay primitive (ODELIA-6), and the
  scalar decision (PROTO-1 → PLANT-1) can proceed in parallel; they converge at
  PLANT-3/4a.
- odelia-native tests (ODELIA-4/4a/4b) track their features and gate the odelia
  release independently of plant; plant's oracle (UX-2) gates the plant ports.
- FF16 invasion (PLANT-4) is the first end-to-end proof and should land before the
  resident path (PLANT-4a) and any TF24 work.
- Nothing merges without UX-2 (the AD-vs-AD oracle) green.

## Not in this release
- Second-order / Hessian (`fwd_adj`).
- **Adaptive sub-stepping in the replay** — the future hardening that would extend the
  TF24f resident coupled gradient past its stiff long-horizon limit (§7.5). v1 gates it
  with a clear error driven by the double replay's env error, never a wrong number.
- Full resident-feedback TF24f at long patch lifetime (the stiffness limit above).
- Moving plant's leaf-level forward-mode AD into odelia (stays plant-local).
