# AD infrastructure: the R boundary for plant and odelia

Companion to [`ad-infrastructure-design.md`](./ad-infrastructure-design.md) and
[`ad-issues.md`](./ad-issues.md). This document covers how AD crosses the R/C++
boundary — the part Rcpp and RcppR6 make awkward — and the interface design that
keeps it clean.

---

## 1. The problem

Rcpp's `as<>`/`wrap<>` and RcppR6's generated bindings marshal **`double`** (and
containers of it). An XAD active type (`xad::adj<double>::active_type`) is a
different scalar with a tape identity; it has **no `as`/`wrap`** and must never
acquire one — serialising it would silently drop the derivative. Three concrete
frictions follow:

1. **Distinct pointer types, hand-dispatched.** `Solver<System<double>>` and
   `Solver<System<active>>` are unrelated C++ types, so they are two different
   `Rcpp::XPtr`s. odelia dispatches between them with a `bool active` flag
   threaded through *every* export (`if (active) get_solver<ActiveSystemType>` …;
   `Solver_reset_impl<SystemType, ActiveSystemType>`). The flag and the pointer's
   real type are kept in sync **only by convention** — `Solver_fit_impl` reads the
   handle as `<ActiveSystemType>` unconditionally, so a double handle reaches it as
   a type-confused reinterpret (UB / crash), not an error.
2. **Two objects leak into R.** To *fit*, the odelia user must construct a
   *second*, active solver (`LeafSolver_new(active = true)`) alongside the double
   one used to simulate. The active/double split becomes the R user's problem.
3. **The plant round-trip.** RcppR6 marshals a `FF16_Environment` to/from an R
   list. The spike's early replay rebuilt each per-stage environment in C++ with
   `Rcpp::as<plant::FF16_Environment>(st[k])` — a serialisation round-trip that is
   lossy for crown-sampled light and rebuilds the whole object per access
   (O(stand), ~1600× slower).

---

## 2. The invariant

> **Only `double` crosses the R boundary. Active types are created, used, and
> destroyed entirely within a single C++ call. R never holds an active handle.**

Everything below follows from this one rule. The AD lifecycle (build active
system → seed inputs → record → run → back-propagate → read adjoints) is an
internal detail of one C++ function; its inputs are doubles (and object handles)
and its outputs are doubles. The only place an active value is read is
`xad::value(x)` / `xad::derivative(x)` → `double`, at the very end.

This is not a limitation to work around — it is the design. It removes the pointer
duplication, the `active` flag, and the round-trip all at once, because none of
them need to exist if active types never travel.

---

## 3. odelia R interface

### 3.1 Retire the R-visible active solver and the `active` flag

R holds only the **double** `Solver`. Gradients are a driver call that builds the
active replay internally:

```cpp
// [[Rcpp::export]]  — R passes the DOUBLE solver handle; only doubles return.
Rcpp::List Solver_gradient(SEXP solver_xp,
                           Rcpp::Nullable<Rcpp::NumericVector> ic,
                           Rcpp::Nullable<Rcpp::NumericVector> params) {
  auto solver = get_solver<SystemType>(solver_xp);          // the DOUBLE solver
  auto [value, grad] = ode::compute_gradient(*solver, independents(ic, params),
                                             sum_of_squares_functional);
  return Rcpp::List::create(_["value"] = value, _["gradient"] = wrap(grad));
}
```

`compute_gradient`/`compute_jacobian` own the active system for the duration of
the call. The `bool active` parameter, the `ActiveSystemType` XPtr, and the dual
`_impl<SystemType, ActiveSystemType>` signatures disappear from the R surface. The
type-confusion hazard is gone because there is only one handle type.

### 3.2 A `rebind` contract so the driver can lift double → active

Today `LeafSolver_new` constructs the active system by hand
(`ActiveSystemType sys_active(pars, drv); sys_active.set_initial_state(…)`). Give
the System a single lift so the driver does this generically:

```cpp
struct System {
  using value_type = S;
  template <class S2> using rebind = System<…, S2>;   // double -> active mould
  template <class S2> System<…,S2> rebind_from() const; // copy config into S2
};
```

The driver calls `rebind_from<active>()`, seeds via the existing `set_params` /
`set_initial_state`, runs, and reads adjoints — all internal.

### 3.3 Tape reuse without exposing active types

The persistent tape matters *within* a call (a Jacobian records once and sweeps
per row) and across calls (an optimizer loop calling `Solver_gradient` repeatedly).
Keep it internal: cache the tape (and optionally the active scratch system) as an
**opaque handle owned by the double Solver**, keyed to nothing R can see. R reuses
the tape by reusing its double Solver; it never learns the tape exists.

### 3.4 No `wrap`/`as` for active types — by policy

Do **not** add an `Rcpp::wrap`/`Rcpp::as` specialization for XAD active types.
Forcing an explicit `xad::value()` / `xad::derivative()` at the one extraction
point keeps derivative loss visible and impossible to do by accident.

---

## 4. plant R interface

### 4.1 Pass the RcppR6 handle; unwrap the pointer, don't serialise

plant's SCM/Patch are RcppR6 objects over the **double** `<T,E>` types; keep them
that way. Gradient entries are hand-written `[[Rcpp::export]]` functions that take
the RcppR6 SCM and unwrap only its pointer — the pattern the spike's best entry
already uses:

```cpp
// [[Rcpp::export]]
Rcpp::NumericMatrix stand_gradient_cpp(SEXP scm_sexp, /* metric names, traits, … */) {
  // Rcpp::as<RcppR6<…>> unwraps .ptr only — NO serialisation of the SCM.
  auto scm = Rcpp::as<plant::RcppR6::RcppR6<plant::SCM<FF16, FF16_Environment>>>(scm_sexp);
  // build the active replay from the LIVE patch (native pointers into
  // environment_history / step_history) — no Rcpp::as<FF16_Environment>.
  // call odelia::compute_jacobian with an EmergentFunctional; return doubles.
}
```

### 4.2 Native harvest — delete the `Rcpp::as<Environment>` round-trip

The active replay reads the resident schedule and per-stage environments **by
pointer into the live Patch's own storage** (issue PLANT-4/5), so the round-trip of
§1.3 never happens. The environment is only ever read as a `double` value plus an
analytic/AD contribution; the active type never touches an RcppR6 object.

### 4.3 One thin R wrapper over one C++ entry

`stand_gradient()` keeps its signature (`scm, metrics, traits, species, feedback`)
but stops branching in R across `native` / `impl` / resident variants. It forwards
to a single C++ entry that dispatches strategy and feedback **in C++** (where the
active types live), returning the double Jacobian + values. RcppR6 never sees an
active type; the emergent translation unit is the only place they exist.

### 4.4 Functionals are constructed in C++, not passed from R

R cannot pass a C++ functional. R passes the **metric names** (strings); plant's
C++ entry maps them to an `EmergentFunctional` (a compile-time object reusing the
scalar-templated reductions) and hands that to odelia's driver. The "functional
shape" is a C++ concept that never crosses to R — only its *selection* does.

---

## 5. Resulting UX

- **odelia:** one solver object; `solver$gradient(params = …)` returns
  `list(value, gradient)`. No `active =` flag, no second solver, no way to hold a
  handle of the wrong type.
- **plant:** `stand_gradient(scm, metrics, traits, species, feedback)` on the
  ordinary (double) SCM the user already ran; returns a double Jacobian. No AD
  objects, no cache marshalling, no round-trip. The replay levels (design §7) are
  internal: the user runs the resident once with the caching their gradient needs,
  and the call selects L0–L3 accordingly.

### 5.1 What the user must run (replay levels, design §7)

The levels are internal, but they determine the *one* thing the user sets on the
resident run:

| Priority | Gradient | User runs | Levels | Canopy |
|---|---|---|---|---|
| primary | census resident (LAI, biomass, basal area) | `control(save_RK45_cache = TRUE)` | L0·L1·L2·L3 | reconstructed (active) |
| next | offspring / invasion | `control(save_RK45_cache = TRUE)` | L0·L1·L3 | frozen |
| advanced | ODE calibration | adaptive run, then fit — no cache; needs observations + loss | L1 | n/a |

`save_RK45_cache = TRUE` is the single AD-relevant control for the emergent
workflows. It records what the replay needs (design §7.5): the resolved ODE step
times, and the adaptive light-spline knot positions / quadrature nodes. The invasion
gradient reuses the recorded resident light frozen; the resident gradient re-runs
`compute_environment` on the recorded knots with active cohorts (no stand-state
cache). A gradient call validates the recording is present and errors clearly if not
(§6.7) — it never silently returns a wrong number. (If the flag name should read as
"prepare for gradients" rather than an implementation detail, that is a small rename
to settle during RIF-7.)

The calibration row is deliberately last: it additionally requires the user to
supply observations and a likelihood, which the emergent workflows do not.
Calibrating the *plant SCM* to data is not a different replay — it is the **resident
replay (L0·L1·L2) plus a likelihood functional**, an addition over the emergent
resident gradient, not a subtraction. The `L1`-only row above is the bare-ODE
(odelia Lorenz) degenerate case; its `set_target`/`advance_target` fit shape
(`test-ad-workflow.R`) serves that case and must **not** set the shape of the primary
`stand_gradient()` UX.

---

## 6. User stories

Each story is the R experience of one persona. The recurring point is what they
*never* touch — the invariant of §2 paying off. The last line of each traces to the
design and surfaces any requirement. Stories are ordered by priority: the emergent
gradients (6.1 resident, 6.2 invasion) are the **primary** plant workflows and need
no observations; calibration (6.3) is an **advanced** case that additionally assumes
targets and a likelihood, and is included to test the invariant under a hot loop —
not as the entry point.

### 6.1 Trait sensitivity — forest ecologist (plant, resident/total)

*"I have a calibrated FF16 stand; how do emergent LAI and biomass respond to leaf
mass per area and wood density?"*

```r
scm <- run_scm(params, control = control(save_RK45_cache = TRUE))
g   <- stand_gradient(scm, metrics = c("LAI", "biomass"),
                      traits = c("lma", "rho"), feedback = "resident")
g$jacobian     # 2x2 doubles: d(metric)/d(trait), with self-shading feedback
```

Never touches XAD, an active type, a tape, or a second solver. The only AD-aware
line is `save_RK45_cache = TRUE`, which belongs to the model run. *Traces to §4,
§5; requires the native harvest (RIF-6) and C++ functional-by-name (§4.4).*

### 6.2 Selection gradient — evolutionary ecologist (plant, mutant/invasion)

*"Give me the invasion-fitness gradient of a rare mutant against an established
resident, so I can locate singular strategies."*

```r
resident <- run_scm(resident_params, control = control(save_RK45_cache = TRUE))
sel <- offspring_production_gradient(resident, traits = c("lma", "hmat"))
# named double vector: d(fitness)/d(trait), resident canopy held frozen
```

The frozen canopy and the `run_mutant` replay are entirely under the hood; the user
picks the workflow by choosing the function, not by managing state. *Traces to the
two-workflow model (design §6.1) and §4.3 (C++ dispatch).*

### 6.3 Gradient-based calibration — modeller fitting data (the hot loop) — *advanced*

*"Run L-BFGS over traits to fit observations; call me for value and gradient each
iteration."* (Advanced: unlike 6.1/6.2 this assumes the user has defined
observations and a likelihood.)

```r
solver <- Solver$new(system, control)     # ONE ordinary (double) solver
solver$set_target(times, obs, obs_idx)
optim(par,
      fn = \(p) solver$value_and_gradient(p)$value,     # doubles in/out
      gr = \(p) solver$value_and_gradient(p)$gradient,
      method = "L-BFGS-B")
```

This is the **L1** replay (design §7): an adaptive run captures the ODE step
schedule, then AD replays pinned to it — no environment cache. It is the story that
most stresses the invariant — a tight loop that would be the natural place to "cache
the active solver in R." It does not: the tape is reused through the opaque cache on
the double solver (§3.3), and the user never sees an active handle or an `active =`
flag. **Requirement surfaced:** expose a single `value_and_gradient(p)` that returns
both from **one** recording, so `fn`/`gr` don't each re-run the tape (a pure
`gradient()` and `loss()` would double the work in an optimizer). *Traces to
§3.1/§3.3, RIF-1/RIF-3; adds `value_and_gradient` to RIF-1.*

### 6.4 A new emergent metric — plant developer

The shipped metrics are **LAI, biomass, and basal area** (plus
`offspring_production`). Extensibility is the story: a developer adds "mean canopy
height" (a tutorial-level example) with one scalar-templated kernel that reuses
existing model functions and registers its name. Then:

```r
stand_gradient(scm, metrics = "mean_height", traits = ff16_default_traits())
```

works immediately — the C++ entry maps the name to the functional; the reverse
sweep is odelia's. No tape code, no odelia change. *Traces to §4.4 and design §5;
confirms "adding a metric is one kernel." mean-height lives in a tutorial, not the
shipped set.*

### 6.5 A new ODE model — odelia downstream developer

*"I'm building a hydraulics module on odelia; I want gradients without writing
marshalling or tape code."*

The developer implements the System contract — `value_type`, `set_params`,
`set_initial_state`, and `rebind` (§3.2) — and gets `Solver_gradient` /
`Solver_jacobian` for free, doubles only. odelia never learns what a hydraulic
segment is. *Traces to §3.2; validates the plant-agnostic odelia surface.*

### 6.6 Trusting a gradient — maintainer (verification)

*"Check the AD gradient against a finite difference before I rely on it."*

```r
g_ad <- stand_gradient(scm, "offspring_production", "lma")$jacobian
g_fd <- (op_at(lma + h) - op_at(lma - h)) / (2 * h)   # re-run scm, perturbed lma
stopifnot(abs(g_ad - g_fd) < tol)
```

Because the AD path returns doubles in the same shape as the finite-difference path,
verification is a plain numeric comparison. *Traces to the doubles-only boundary and
the regression oracle (UX-2).*

### 6.7 Edge case — forgot the cache (fail loud, never confuse)

*"I ran the SCM without `save_RK45_cache` and asked for a gradient."*

```r
scm <- run_scm(params)                 # no cache
stand_gradient(scm, "LAI", "lma")
#> Error: no resident schedule cached; re-run with control(save_RK45_cache = TRUE)
```

A clear, actionable error — never a crash. This is the deliberate contrast with the
§1.1 hazard: with active types off the R surface there is no wrong-typed handle to
reinterpret, so the boundary can only fail loudly. *Traces to §2, §3.1.*

---

## 7. Work items (feed `ad-issues.md`)

| Item | Repo | Class | Note |
|---|---|---|---|
| RIF-1 | odelia | CP | `Solver_gradient`/`Solver_jacobian`/`value_and_gradient` on the double handle; retire the `active` flag + `ActiveSystemType` XPtr from the R surface |
| RIF-2 | odelia | CP | `rebind` lift contract (§3.2) so the driver constructs the active system |
| RIF-3 | odelia | CP-support | opaque tape/active-scratch cache on the double Solver (§3.3) |
| RIF-4 | odelia | NTH | policy check: ensure no `wrap`/`as` for active types compiles |
| RIF-5 | plant | CP | single `stand_gradient_cpp` entry taking the RcppR6 handle; C++ strategy/feedback dispatch (§4.1, §4.3) |
| RIF-6 | plant | CP | native-pointer harvest; delete `Rcpp::as<*_Environment>` (ties to PLANT-4/5) |
| RIF-7 | plant | CP-support | thin `stand_gradient()` R wrapper; remove R-side branching |

**Dependencies:** RIF-1..3 sit under ODELIA-1/2/3 (they are the R-facing side of the
driver). RIF-5..7 sit under PLANT-4/5 (native harvest) and UX-1 (stable surface).
None require exposing an active type to R; that is the point.

## Open questions

- **Keep any R-visible active solver at all?** An advanced user driving a bespoke
  optimizer loop in R might want to hold AD state across calls. The opaque-cache
  approach (§3.3) covers the common case; decide whether a power-user escape hatch
  is worth the re-introduced footgun. Recommendation: no — start with the invariant,
  add an escape hatch only on a demonstrated need.
- **RcppR6 template-type dispatch.** `stand_gradient_cpp` must resolve the SCM's
  `<T,E>` (FF16 / TF24 / TF24f) to instantiate the right replay. Confirm whether
  this dispatch is best done by the existing `extract_RcppR6_template_types` (R
  side, strings) feeding a C++ `switch`, or a C++-side type tag on the handle.
