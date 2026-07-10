# AD infrastructure design: plant emergent gradients on odelia's AD runtime

**Scope:** `traitecoevo/plant` (SCM emergent-trait gradients) on `traitecoevo/odelia`
(the AD-aware ODE runtime). A validated prototype — plant PR
[#553](https://github.com/traitecoevo/plant/pull/553), tracking issue
[#472](https://github.com/traitecoevo/plant/issues/472) — is the specification and the
regression oracle.
**Companions:** [`ad-record-replay.md`](./ad-record-replay.md) — the record→replay
primitive beneath the replay levels (§7); [`ad-r-interface.md`](./ad-r-interface.md) —
the R/C++ AD boundary for both packages.

---

## 1. Thesis

odelia is the AD runtime for this package family: it compiles XAD's `Tape` once and
ships a scalar-templated ODE `Solver`, a persistent tape, a reverse-mode gradient
driver, and a differentiable spline. plant already `LinkingTo: odelia` and runs its SCM
on `odelia::ode::Solver<patch_type>`.

Exact reverse-mode SCM gradients therefore follow from **surgical changes to the plant
components that already exist**, plus a **small generic surface on odelia** — not a
second, plant-private AD stack. Two facts carry the design:

- **The Patch is already the odelia System** (`Solver<patch_type>`); it needs only to
  satisfy odelia's AD contract.
- **`SCM::run_mutant()` is already the frozen-schedule replay** — it pins integration to
  the resident's cached `step_history`, reads the cached environment, and suppresses the
  mutant's self-competition. The invasion gradient is the derivative *of that existing
  run*, not of a new engine.

---

## 2. Decisions

1. **Reverse mode at the SCM level.** The emergent Jacobian is metrics × traits
   (≈4 × 28) — outputs ≪ inputs — so the reverse (adjoint) driver is optimal.
2. **One uniform `value_type = S`.** Patch/Species/Node/Individual run a single active
   scalar throughout; there is no separate active/frozen physiology axis. The scalar
   lives with the types that own the parameters, so `<T,E,S>` collapses to the existing
   `<T,E>` (§5.1).
3. **Surgical changes to existing types; no new abstractions.** No `SCMSystem`, no
   `StrategyConcept`. Strategy and Patch already exist and are made AD-compatible in
   place (§5).
4. **odelia stays plant-agnostic.** It knows Systems, tapes, functionals, and Jacobians
   — never cohorts, traits, emergence, or basal area. plant leads odelia's API by
   concrete need.
5. **Two workflows, both differentiated.** *Resident/total* (the stand differentiates
   with self-feedback) and *mutant/invasion* (a rare mutant in the frozen resident
   canopy) already exist as `run()` and `run_mutant()`; the design differentiates each
   (§6).
6. **Forward mode stays plant-local.** The leaf gas-exchange optimizer keeps its
   `xad::fwd` IFT solve; odelia accommodates it through a supplied-derivative edge
   (§4, §5.2), not by absorbing it.

---

## 3. Three layers

The work separates cleanly, and the document is organised by these layers.

| Layer | Owns | Knows about plant? |
|---|---|---|
| **odelia** (§4) | tape, Solver, `compute_gradient`/`compute_jacobian`, the functional *shape*, the differentiation-target *shape* | No |
| **plant** (§5) | Patch-as-System, Strategy AD-compat, scalar-templated reductions, differentiating `run`/`run_mutant`, SCM orchestration | — |
| **UX / API / workflow** (§6) | `stand_gradient()` and friends, the two workflows, the metric set | — |

---

## 4. odelia layer (generic, plant-agnostic)

odelia supplies the compiled single-definition `Tape`, a scalar-templated
`Solver<System>` on `System::value_type` with a Solver-owned persistent tape, the
`set_params`/`set_initial_state` System contract with the `rebind` lift and `scatter`
routing, and the differentiable spline — over the vendored XAD facilities
(`computeJacobian`, `CheckpointCallback`, `xad::adj`/`xad::fwd`). The generic AD surface
above them is four pieces.

**(a) The functional shape.** A functional maps a replayed System to output scalar(s),
and nothing more. It is a **pure reduction**: it reads state and returns a scalar; it
does not drive the solver and carries no schedule. The driver owns the replay on the
recorded schedule (`advance_fixed(recorded_steps())`) and the Solver stores no fit
state. The calibration loss is one prebuilt instance, `least_squares`, holding only its
measured observations and which recorded steps they attach to; plant supplies its own
emergent functional (§6). odelia never learns what the scalars mean.

```cpp
// F is any callable  std::vector<S>  f(const System& replayed).
template <class System, class F>
std::pair<double, std::vector<double>>
compute_gradient(Solver<System>&, const DifferentiationTargets&, F&& functional);
```

**(b) `compute_jacobian` (reverse, generic).** Record once, one adjoint sweep per output
row — the adjoint variant of `XAD::computeJacobian`, over any odelia System, with no
plant knowledge. Two properties are load-bearing:

- **Pass the codomain** (output count) to `computeJacobian`. Omitting it costs an extra
  forward evaluation — a full solve, i.e. a full SCM replay for plant — just to size the
  outputs. The functional knows its output arity, so it supplies it.
- **Column order is a contract.** Output column `j` is `d(output)/d(leaf_j)` for the
  `j`-th seeded leaf *in the order the System consumes them*. A caller resolving trait
  names → fields must pack and scatter in the same order, or the columns transpose
  silently.

**(c) `DifferentiationTargets` — the "gradients w.r.t. what?" shape.** A generic
description of which registered inputs are active, with no plant vocabulary: a set of
addressed leaves the System routes to its own fields. A trait-seed and an
initial-condition-seed are the same kind of thing — a registered leaf — and IC
sensitivity is a genuine target (plant's `make_initial_state` / `export_patch_state`
make an initial size distribution a first-class input). So the shape carries no
privileged `params` / `initial_state` split: any subset — traits, ICs, or both — is
expressible, and the System owns the routing. odelia never learns what a leaf means.

**(d) `SuppliedDerivative` — accommodating IFT / forward-mode results.** For a value the
forward pass computes off-tape (a root-find or optimizer result), its partials are known
analytically but were never recorded. `SuppliedDerivative` injects them into the reverse
tape via `CheckpointCallback`: it registers the off-tape value as a fresh leaf and, on
the reverse sweep, increments each input's adjoint by `ybar · ∂y/∂x_i`. This is odelia's
compatibility seam for plant's leaf-optimizer forward-mode sensitivity and the TF24
stomatal IFT — generic (odelia sees inputs, an output, and partials). It is a free
function called from within the forward pass, where the off-tape value actually exists,
not a pre-declared target.

---

## 5. plant layer (surgical changes to existing components)

Each change below is a small in-place modification of a type that already exists.

### 5.1 Patch — satisfy odelia's System AD contract

The Patch is already the System (`SCM` holds `odelia::ode::Solver<patch_type>`).

- **Uniform `value_type = S`.** Let the strategy and environment carry the scalar
  (`T = FF16_Strategy<S>`, `E = FF16_Environment<S>`) and have
  `Individual`/`Node`/`Species`/`Patch` derive `using value_type = typename
  T::value_type`. That keeps the existing `<T,E>` template shape, matches odelia's
  `value_type` convention, and puts the scalar with the types that own the parameters.
- **`set_params(tape, it)` / `set_initial_state(tape, it, t0)`** returning the registered
  active inputs — the same contract `leaf_thermal_system.hpp` models. This is where
  traits and birth-rate become active; the map from trait names to registered fields is
  plant's (§5.2).

These are methods on the existing Patch, not a new wrapper type.

### 5.2 Strategy — make traits seedable, keep physiology scalar-templated

The Strategy classes already expose per-method `template<class S>` physiology
(`area_leaf<S>`, `update_dependent_aux<S>`). Two changes:

- **One parameter representation that can be active.** The strategy stores a single
  scalar-templated parameter store, so a named trait registers active directly (feeding
  §5.1's `set_params`) — no separate double `pars` plus a lifted active struct.
- **Leaf-optimizer edge.** TF24/TF24f keep their forward-mode leaf solve; its sensitivity
  crosses into the reverse sweep as a `SuppliedDerivative` (§4d).

FF16 (no optimizer) needs only the parameter-store change.

### 5.3 Reductions — scalar-template the model's own quadratures

`Species::compute_competition(height)` **is** the census/light trapezium. Template it and
the patch census integral on `S`, so one function serves the forward model (`S=double`,
unchanged) and the gradient replay (`S=active`). Exactness is then structural — the
gradient reads the model's own reduction, not a maintained bit-for-bit copy.

### 5.4 The replay is `run` / `run_mutant`, differentiated

The two workflows differentiate two existing runs, differing precisely in how the canopy
(L3, §7) is accessed:

- **Invasion** = differentiate `SCM::run_mutant(p)` with `S=active`. It switches the
  Patch to the cached resident environment (`set_mutant()`), pins the schedule to the
  resident `step_history`, and runs — the mutant reads the **frozen** canopy and does not
  compete with itself (`is_mutant_run` gates `compute_competition`). The derivative
  through the canopy is zero.
- **Resident / total** = differentiate the resident run on the *same* frozen L0/L1
  schedule, but with the canopy **recomputed live from the active cohorts** so a trait
  re-shades the stand. It must **not** read the frozen `environment_history` — that
  collapses it to the invasion gradient.

The resident canopy is re-computed by re-running `compute_environment` on the recorded
fixed knots (§7), so the SCM is differentiated directly and no per-RK-stage stand-state
cache is kept.

### 5.5 SCM orchestration

The SCM stays the plant-specific orchestrator; odelia supplies the atoms (`advance_fixed`
/ `advance_adaptive`, tape). Two orchestration points matter for AD:

- **Node introductions grow the system; they are not discontinuities.** When the node
  schedule is frozen (L0), the introduction times are *constants*
  (`d(t_intro)/d(trait) = 0`), so introductions add tape variables (a wider state) but
  inject no discontinuity into the differentiated output. This is exactly why the
  schedule must be frozen: an adaptive, trait-dependent schedule would be
  non-differentiable. One recording spans the whole run across introductions.
- **Adaptivity is frozen in layers, not all at once.** The SCM has four adaptive
  constructions that each break differentiability, frozen independently at different
  depths — the replay levels of §7.

---

## 6. UX / API / workflow layer

### 6.1 The two workflows

| Workflow | SCM entry | Environment | Gradient meaning |
|---|---|---|---|
| **Resident / total** | `run()` | co-moving (stand re-shades itself) | d(emergent metric)/d(trait), full self-feedback |
| **Mutant / invasion** | `run_mutant()` | frozen resident canopy | selection gradient: d(rare-mutant fitness)/d(mutant trait) |

The mutant has low density in an established stand: it competes with the residents but
not with itself, so its fitness gradient holds the environment fixed (positive fitness ⇒
it can invade). These are two genuinely separate workflows — pre-existing SCM
capabilities the design differentiates — not a mode flag on one.

### 6.2 Gradient semantics

| Gradient | Workflow | Cross term | Notes |
|---|---|---|---|
| Trait, invasion | mutant | 0 (frozen canopy) | the selection gradient; `offspring_production` is always this |
| Trait, resident | resident | present | species re-shades the stand it lives in |
| Birth-rate, census metric | resident | present | the frozen part is the identity `metric / birth_rate`; the resident (canopy-feedback) axis needs a tape and is non-trivial — it **flips the sign of biomass** vs. the identity |
| Birth-rate, `d R0/d birth_rate` | mutant framing | — | `d(net_reproduction_ratio)/d(birth_rate)`; the density-feedback axis |

**The demographic-equilibrium solve.** `d R0/d birth_rate` is the plant-side derivative
for the Newton solve that finds the equilibrium birth rate (R0 = 1) — the resident
density at which a strategy sustains itself. It is computed on the same coupled replay via
the mutant framing (mutant traits = resident; the change in mutant fitness as resident
density moves is the density feedback). The birth-rate gradient is therefore not only a
sensitivity but the enabling piece for equilibrium and invasion analyses.

So the resident-vs-frozen distinction is never a minor correction: the canopy feedback can
dominate and reverse the sign of a response. "Which feedback" is a modelling choice, not a
tolerance detail.

### 6.3 Public API

`stand_gradient(scm, metrics, traits, species, feedback)` keeps its shape and maps to
`compute_jacobian(solver, targets(traits, species, birth_rate),
EmergentFunctional{metrics})`, where `EmergentFunctional` is plant-supplied and reuses
§5.3's reductions. The R harvest and `Rcpp::as<>` round-trip disappear — the functional
reads native state on the tape. Adding a metric is a one-kernel change in plant that
reuses model functions; odelia is untouched.

### 6.4 Two orthogonal axes

A workflow is a choice on **two independent axes**, and separating them removes most of
the apparent complexity:

- **Replay** — which adaptive constructions must be recorded and replayed fixed (§7). A
  property of the *system*, not the question: the plant SCM always needs L0·L1·L2; a bare
  ODE needs only L1. Plus the ecological choice of **feedback** — resident (canopy
  re-run, active) vs. invasion (canopy frozen), which is L3.
- **Functional** — what scalar is differentiated: an emergent metric (no observations) or
  a likelihood over observations. Orthogonal to replay.

| Priority | Workflow | Replay (system) | Feedback | Functional |
|---|---|---|---|---|
| **1 (primary)** | Emergent gradient, resident | SCM: L0·L1·L2 | resident (re-run) | emergent metric |
| **2** | Mutant / invasion fitness | SCM: L0·L1·L2 | invasion (frozen) | emergent metric |
| **3 (advanced)** | Calibration / inference | same replay as the system needs | resident | likelihood over observations |

**Calibration is a different functional, not a different replay.** Calibrating the plant
SCM to data uses the *same* resident replay as workflow 1 with a likelihood functional on
top — an addition, not a subtraction. The bare-ODE fit (odelia Lorenz, L1 only) is the
degenerate case; its shape must not set the plant UX. The emergent workflows are primary
because they need no observations; calibration is advanced because it additionally
requires the user to define observations and a likelihood.

---

## 7. Replay levels: four independent freezes

The one record→replay idea (mechanism in
[`ad-record-replay.md`](./ad-record-replay.md)) applies at four *distinct* depths in the
SCM. Each freezes a different adaptive construction that would otherwise break
differentiability; they compose, and a given gradient needs only some.

| Level | Freezes | Removes non-diff from | Owner |
|---|---|---|---|
| **L0 — node schedule** | which cohorts exist and when introduced | adaptive cohort introduction | plant |
| **L1 — ODE step times** | the adaptive RKCK step selection (replay `advance_fixed`) | adaptive step-size control | **odelia** |
| **L2 — quadrature / interpolator knots** | height-quadrature abscissae and light-spline knots | adaptive quadrature / interpolation refinement | plant + odelia spline |
| **L3 — resident canopy** | the resident environment the focal cohorts read | canopy feedback (frozen for a mutant) | plant |

L1 alone is the bare-ODE calibration case: run adaptively, capture `times()`, replay
`advance_fixed(recorded_steps())` — no environment cache, just the resolved schedule.

**L2 has two variants.** For the resident light spline, the knots are frozen (positions)
and the values active — odelia's differentiable spline. For a census integrated over
height the integration bound *is* an active plant height, so the quadrature **nodes
move**: this needs the scalar-templated `QK`, not a frozen-node replay, which would miss
the moving-node sensitivity.

**L3 is the correctness crux — frozen ≠ reconstructed.** The resident run records the
resident environment *values* and the light-spline knot *positions*; the two gradients
diverge on which they read.

- **Invasion (frozen):** a rare mutant reads the recorded environment as `double`
  background — off the tape entirely, so the derivative through the canopy is zero by
  construction (not an active constant with a zeroed derivative). This is `run_mutant`.
- **Resident / total (re-run):** the canopy is recomputed live by re-running
  `compute_environment` on the recorded knots with the active, re-evolved cohorts, so a
  trait re-shades the stand through `area_leaf`. **Replaying the frozen environment here
  would silently give the invasion gradient** — the self-shading cross term would be
  missing. Only the knot *positions* are recorded; the cohort values come from the replay
  itself, so no per-RK-stage stand-state cache is needed.

### Which gradient needs which levels

| Gradient | L0 | L1 | L2 | L3 |
|---|:--:|:--:|:--:|:--:|
| Bare-ODE fit (odelia Lorenz — degenerate) | | ✅ | | |
| Offspring, invasion | ✅ | ✅ | | frozen |
| Census (LAI/biomass/basal area), resident | ✅ | ✅ | ✅ | **re-run** |
| Census, invasion | ✅ | ✅ | ✅ | frozen |
| SCM calibration to data | ✅ | ✅ | ✅ | re-run |

The last two rows share the L0·L1·L2·L3-re-run replay and differ only in the functional
(emergent metric vs. likelihood) — the point of §6.4.

**The fixed comb keeps functionals simple.** With L0–L1 frozen, the emergent functionals
are reductions over a fixed set of cohorts. Offspring is a constant-weighted sum
`Σ tw_i · offspring_i` (fixed introduction times ⇒ constant weights). A census reduction
is `Species::compute_competition` scalar-templated (§5.3), with L2 supplying the (moving
or frozen) quadrature nodes. The design differentiates the *same* quadratures the model
already uses; it adds no new integration.

---

## 8. Strategy coverage and the cross-sensitivity term

| Strategy | Offspring | Census, invasion | Census, resident |
|---|---|---|---|
| FF16 (no leaf optimizer) | ✅ | ✅ | ✅ |
| TF24 / TF24f (leaf optimizer) | ✅ | ✅ | ✅ (via the cross-term below) |

TF24/TF24f solve a per-cohort leaf optimum, and a census metric depends on it: a trait
shift moves the optimum → growth → census density. The design captures this because the
whole `run_mutant` replay runs on **one tape** with the leaf IFT delivered as a
`SuppliedDerivative` (§4d): the reverse sweep traverses density→optimum→trait natively,
rather than a linearised, frozen harvest that injects the leaf sensitivity along the
focal path but zeroes the cross-term through the density. This full-scalar path is what
puts TF24 census on the same footing as FF16.

---

## 9. Boundaries

- **Stiff TF24f resident coupling at long patch lifetimes.** Fixed-node replay assumes
  the recorded nodes stay adequate when the scalar goes active. This holds for FF16 and
  for TF24/TF24f at moderate horizons. It fails for the TF24f resident coupled feedback at
  long lifetimes: the `log_density`↔canopy loop is stiff, and the live SCM tames it with
  *adaptive* sub-stepping that a fixed schedule cannot reproduce — the frozen-step replay
  itself drifts. The true gradient exists (full-SCM finite differences are finite); only
  the fixed replay diverges. The primitive is the right default; the stiff coupled horizon
  needs adaptive sub-stepping in the replay, and until it has that the case is gated by a
  clear error (driven by the double replay's environment error), never a wrong number.
- **Zero-height cohort trap.** A cohort introduced on the final step (`birth == N`) can
  sit at zero height, where `area_leaf = (h/a_l1)^(1/a_l2)` differentiates to
  `0·log(0) = NaN` (it also biases the value). Establishing `birth ≥ N` cohorts at the
  seed height `h0` fixes both — a concrete AD trap that carries an explicit test.
- **Second order (Hessian) is out of scope.**
