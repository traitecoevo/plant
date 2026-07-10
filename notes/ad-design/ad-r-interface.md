# AD infrastructure: the R boundary for plant and odelia

Companion to [`ad-infrastructure-design.md`](./ad-infrastructure-design.md). This
document covers how AD crosses the R/C++ boundary — the part Rcpp and RcppR6 make
awkward — and the interface design that keeps it clean.

---

## 1. The problem

Rcpp's `as<>`/`wrap<>` and RcppR6's bindings marshal **`double`** and containers of it.
An XAD active type (`xad::adj<double>::active_type`) is a different scalar with a tape
identity; it has **no `as`/`wrap`** and must never acquire one — serialising it would
silently drop the derivative. Three frictions follow, and the invariant of §2 removes all
three at once:

1. **Two pointer types.** `Solver<System<double>>` and `Solver<System<active>>` are
   unrelated C++ types, so a binding that exposes both needs a flag to dispatch between
   them — and the flag and the pointer's real type stay in sync only by convention, a
   type-confusion footgun.
2. **Two objects in R.** Exposing the active solver forces the user to build and hold a
   *second* solver alongside the double one they simulate with.
3. **The plant round-trip.** Marshalling a `FF16_Environment` to/from an R list to
   rebuild each per-stage environment is lossy for crown-sampled light and rebuilds the
   whole object per access — O(stand), orders of magnitude slower than a native-pointer
   read.

---

## 2. The invariant

> **Only `double` crosses the R boundary. Active types are created, used, and destroyed
> entirely within a single C++ call. R never holds an active handle.**

Everything below follows from this one rule. The AD lifecycle — build the active system,
seed inputs, record, run, back-propagate, read adjoints — is an internal detail of one
C++ function; its inputs are doubles (and object handles) and its outputs are doubles. The
only place an active value is read is `xad::value(x)` / `xad::derivative(x)` → `double`, at
the very end.

This is the design, not a limitation worked around: the pointer duplication, the `active`
flag, and the round-trip do not need to exist if active types never travel.

---

## 3. odelia R interface

### 3.1 R holds only the double Solver

A gradient is a driver call that builds the active solver internally:

```cpp
// [[Rcpp::export]]  — R passes the DOUBLE solver handle; only doubles return.
Rcpp::List Solver_gradient(SEXP solver_xp,
                           Rcpp::Nullable<Rcpp::NumericVector> ic,
                           Rcpp::Nullable<Rcpp::NumericVector> params) {
  auto solver = get_solver<SystemType>(solver_xp);          // the DOUBLE solver
  auto [value, grad] = ode::compute_gradient(*solver, differentiation_targets(ic, params),
                                             emergent_functional);
  return Rcpp::List::create(_["value"] = value, _["gradient"] = wrap(grad));
}
```

`compute_gradient` / `compute_jacobian` own the active system for the duration of the
call. There is one handle type, so a gradient call cannot reinterpret a double handle as
active — the boundary can only fail loudly.

### 3.2 A `rebind` contract so the driver lifts double → active

The System carries a single lift so the driver builds the active system generically:

```cpp
struct System {
  using value_type = S;
  template <class S2> using rebind = System<…, S2>;     // double -> active mould
  template <class S2> System<…,S2> rebind_from() const; // copy config into S2
};
```

The driver calls `rebind_from<active>()`, seeds via the existing `set_params` /
`set_initial_state`, runs, and reads adjoints — all internal. Only values cross in the
lift, so the active system starts free of tape identity.

### 3.3 Tape reuse without exposing active types

The persistent tape matters *within* a call (a Jacobian records once and sweeps per row)
and *across* calls (an optimizer loop calling `Solver_gradient` repeatedly). It stays
internal: the active solver is cached on the double Solver, and R reuses it by reusing its
double Solver, never learning the tape exists. Three facts:

- **The active solver is the only cached thing.** The gradient runs on it, and a `Solver`
  carries its own `tape`, so that tape *is* the reused tape — there is nothing else to
  cache. It is held on the double Solver, and its type is **named, not erased**: the
  System supplies `rebind`, so `System::rebind<active>` spells it from inside
  `Solver<System>` — no `void*`, no `static_cast`.
- **Anchored on the `Solver` object, not an R handle.** plant's SCM holds the solver as a
  plain C++ member (`scm.h`: `Solver<patch_type> solver;`) and never wraps it in an XPtr,
  so a handle-slot anchor would be invisible to it. Anchoring on the `Solver` object gives
  `stand_gradient_cpp` the reuse for free.
- **The recording is read per call, not frozen into the active solver.** The recording —
  resolved step times and any interpolator/quadrature spacing or frozen field values (the
  replay levels of design §7) — is per-run state owned by the immutable double Solver and
  handed over on every call. Its validity domain is the ICs and params of the double run:
  a replay may vary the mutant / observations / functional, but changing ICs or params
  invalidates the recording and forces a re-record.

Two senses of "cache" stay distinct: the active solver + tape is a **speed** cache; the
record/replay recording is **semantic**. Reusing the cache never changes a number.

### 3.4 No `wrap`/`as` for active types — by policy

Do **not** add an `Rcpp::wrap` / `Rcpp::as` specialization for active types. Forcing an
explicit `xad::value()` / `xad::derivative()` at the one extraction point keeps derivative
loss visible and impossible to do by accident.

---

## 4. plant R interface

### 4.1 Pass the RcppR6 handle; unwrap the pointer, don't serialise

plant's SCM/Patch are RcppR6 objects over the **double** `<T,E>` types; they stay that
way. Gradient entries are hand-written `[[Rcpp::export]]` functions that take the RcppR6
SCM and unwrap only its pointer:

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

### 4.2 Native harvest — no `Rcpp::as<Environment>` round-trip

The active replay reads the resident schedule and per-stage environments **by pointer into
the live Patch's own storage**, so the round-trip of §1 never happens. The environment is
only ever read as a `double` value plus an analytic/AD contribution; the active type never
touches an RcppR6 object.

### 4.3 One thin R wrapper over one C++ entry

`stand_gradient()` keeps its signature (`scm, metrics, traits, species, feedback`) and
forwards to a single C++ entry that dispatches strategy and feedback **in C++**, where the
active types live, returning the double Jacobian and values. RcppR6 never sees an active
type; the emergent translation unit is the only place they exist.

### 4.4 Functionals are constructed in C++, not passed from R

R cannot pass a C++ functional. R passes the **metric names** (strings); plant's C++ entry
maps them to an `EmergentFunctional` (a compile-time object reusing the scalar-templated
reductions) and hands that to odelia's driver. The functional shape never crosses to R —
only its *selection* does.

---

## 5. Resulting UX

- **odelia:** one solver object; `solver$gradient(params = …)` returns `list(value,
  gradient)`. No `active =` flag, no second solver, no way to hold a handle of the wrong
  type.
- **plant:** `stand_gradient(scm, metrics, traits, species, feedback)` on the ordinary
  (double) SCM the user already ran; returns a double Jacobian. No AD objects, no cache
  marshalling, no round-trip. The replay levels are internal.

### 5.1 What the user sets on the resident run

The replay levels are internal, but they determine the one thing the user sets:

| Priority | Gradient | User runs | Levels | Canopy |
|---|---|---|---|---|
| primary | census resident (LAI, biomass, basal area) | `control(save_RK45_cache = TRUE)` | L0·L1·L2·L3 | reconstructed (active) |
| next | offspring / invasion | `control(save_RK45_cache = TRUE)` | L0·L1·L3 | frozen |
| advanced | ODE calibration | adaptive run, then fit (needs observations + likelihood) | L1 | n/a |

`save_RK45_cache = TRUE` is the single AD-relevant control for the emergent workflows. It
records what the replay needs: the resolved ODE step times and the adaptive light-spline
knot positions / quadrature nodes. The invasion gradient reuses the recorded resident
light frozen; the resident gradient re-runs `compute_environment` on the recorded knots
with active cohorts. A gradient call validates the recording is present and errors clearly
if not (§6.7) — it never silently returns a wrong number.

Calibration is last because it additionally requires observations and a likelihood.
Calibrating the plant SCM to data is the resident replay plus a likelihood functional — an
addition over the emergent resident gradient, not a subtraction. The bare-ODE (odelia
Lorenz) L1-only fit is the degenerate case, and its shape must not set the primary
`stand_gradient()` UX.

---

## 6. User stories

Each story is the R experience of one persona; the recurring point is what they *never*
touch — the invariant of §2 paying off. The emergent gradients (6.1 resident, 6.2
invasion) are the **primary** plant workflows and need no observations; calibration (6.3)
is an **advanced** case that additionally assumes observations and a likelihood.

### 6.1 Trait sensitivity — forest ecologist (resident/total)

*"I have a calibrated FF16 stand; how do emergent LAI and biomass respond to leaf mass per
area and wood density?"*

```r
scm <- run_scm(params, control = control(save_RK45_cache = TRUE))
g   <- stand_gradient(scm, metrics = c("LAI", "biomass"),
                      traits = c("lma", "rho"), feedback = "resident")
g$jacobian     # 2x2 doubles: d(metric)/d(trait), with self-shading feedback
```

Never touches XAD, an active type, a tape, or a second solver. The only AD-aware line is
`save_RK45_cache = TRUE`, which belongs to the model run.

### 6.2 Selection gradient — evolutionary ecologist (mutant/invasion)

*"Give me the invasion-fitness gradient of a rare mutant against an established resident,
so I can locate singular strategies."*

```r
resident <- run_scm(resident_params, control = control(save_RK45_cache = TRUE))
sel <- offspring_production_gradient(resident, traits = c("lma", "hmat"))
# named double vector: d(fitness)/d(trait), resident canopy held frozen
```

The frozen canopy and the `run_mutant` replay are entirely under the hood; the user picks
the workflow by choosing the function, not by managing state.

### 6.3 Gradient-based calibration — modeller fitting data (the hot loop) — *advanced*

*"Run L-BFGS over traits to fit observations; call me for value and gradient each
iteration."*

```r
solver <- Solver$new(system, control)     # ONE ordinary (double) solver
solver$set_observations(times, obs, obs_idx)
optim(par,
      fn = \(p) solver$value_and_gradient(p)$value,     # doubles in/out
      gr = \(p) solver$value_and_gradient(p)$gradient,
      method = "L-BFGS-B")
```

This is the **L1** replay: an adaptive run captures the ODE step schedule, then AD replays
pinned to it. It is the story that most stresses the invariant — a tight loop that would
be the natural place to "cache the active solver in R." It does not: the tape is reused
through the opaque cache on the double solver (§3.3), and the user never sees an active
handle. A single `value_and_gradient(p)` returns both from **one** recording, so `fn` and
`gr` share the tape rather than each re-running it.

### 6.4 A new emergent metric — plant developer

The shipped metrics are LAI, biomass, basal area, and `offspring_production`. A developer
adds "mean canopy height" with one scalar-templated kernel that reuses existing model
functions and registers its name:

```r
stand_gradient(scm, metrics = "mean_height", traits = ff16_default_traits())
```

works immediately — the C++ entry maps the name to the functional; the reverse sweep is
odelia's. No tape code, no odelia change.

### 6.5 A new ODE model — odelia downstream developer

*"I'm building a hydraulics module on odelia; I want gradients without writing marshalling
or tape code."*

The developer implements the System contract — `value_type`, `set_params`,
`set_initial_state`, and `rebind` (§3.2) — and gets `Solver_gradient` / `Solver_jacobian`
for free, doubles only. odelia never learns what a hydraulic segment is.

### 6.6 Trusting a gradient — maintainer (verification)

*"Check the AD gradient against a finite difference before I rely on it."*

```r
g_ad <- stand_gradient(scm, "offspring_production", "lma")$jacobian
g_fd <- (op_at(lma + h) - op_at(lma - h)) / (2 * h)   # re-run scm, perturbed lma
stopifnot(abs(g_ad - g_fd) < tol)
```

Because the AD path returns doubles in the same shape as the finite-difference path,
verification is a plain numeric comparison.

### 6.7 Edge case — forgot the cache (fail loud, never confuse)

*"I ran the SCM without `save_RK45_cache` and asked for a gradient."*

```r
scm <- run_scm(params)                 # no cache
stand_gradient(scm, "LAI", "lma")
#> Error: no resident schedule cached; re-run with control(save_RK45_cache = TRUE)
```

A clear, actionable error — never a crash. With active types off the R surface there is no
wrong-typed handle to reinterpret, so the boundary can only fail loudly.

### 6.8 Enabling AD on a system with adaptive numerics — odelia downstream developer

*"My hydraulics module's rates read an adaptively-refined interpolator — a bare ODE like
Lorenz didn't. What must I implement so reverse-mode gradients are correct, and what if I
also want a cheap run that holds part of the system fixed?"*

6.5 answered this for a bare ODE: implement `value_type` / `set_params` /
`set_initial_state` / `rebind`, and gradients come for free. Once the system has adaptive
sub-numerics there are **two further, independent capabilities**, and keeping them separate
is what removes the "cache" / "environment" confusion:

**1. Node-position replay — required for AD, automatic, not a user flag.** Any adaptive
construction makes parameter-dependent discrete decisions about where to place nodes;
differentiating through those branches corrupts the tape. So the double pass records the
node **positions** and the AD pass replays them **fixed**, while the values at those nodes
stay fully differentiable. The ODE step schedule is universal and odelia owns it; if your
system builds an interpolator or quadrature, you record its node positions through the
`Replayable` hooks. This is switched on by *doing reverse-mode AD*, not by any user control
— miss it and gradients are silently wrong wherever the adaptive component bites.

**2. Value freeze — an optional variant workflow.** Separately, you may want a run that
holds some *recomputable quantity* constant — its derivative zero by construction. plant's
rare mutant reading a fixed resident field is one instance, but the quantity need not be an
"environment" or an interpolator; it could be a held-fixed sub-state or frozen state
variables. You record it on the double pass and read it frozen on a designated replay. This
*is* a per-call choice, but the user expresses it by calling the variant entry (a
`run_mutant`-style function), not a `feedback` flag — the workflow function *is* the choice.

odelia stays agnostic to both: it records "some node positions" and "some values" and never
learns a node is a height or a value an environment. You implement the `Replayable` hooks
(`record_stage` / `record_ode_step` / `replay_step` / `has_recorded_field`); the variant
entry (`run` vs `run_mutant`) chooses whether the frozen-field cache is populated, and
`has_recorded_field()` reports it. The *user* does nothing for capability 1 (it rides the
AD call) and picks the variant function for capability 2.

---

## 7. Two design decisions to record

- **No R-visible active solver — not even a power-user escape hatch.** An advanced user
  driving a bespoke optimizer loop in R might want to hold AD state across calls, but the
  opaque cache (§3.3) already covers the hot loop without re-introducing the wrong-typed
  handle. The invariant wins; an escape hatch waits for a demonstrated need.
- **RcppR6 template-type dispatch.** `stand_gradient_cpp` resolves the SCM's `<T,E>`
  (FF16 / TF24 / TF24f) to instantiate the right replay. The dispatch is done in C++ — a
  strategy tag off the handle feeding a `switch` — so the active types never surface in R.
