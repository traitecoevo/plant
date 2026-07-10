# AD record → replay: the one adaptive-numerics primitive

**Scope:** `traitecoevo/odelia` — the record→replay primitive beneath every adaptive
construction on the AD path. It applies in plant unchanged in shape: the light spline
and the crown quadrature record the same way, and the plant `QK` node-set slots into
the identical machinery.
**Companions:** [`ad-infrastructure-design.md`](./ad-infrastructure-design.md) (the
plant emergent-gradient design this primitive sits under);
[`ad-r-interface.md`](./ad-r-interface.md) (the R/C++ boundary and the reuse of the
active solver across calls).

---

## 1. Thesis

Every adaptive numerical construction in the stack — the RKCK stepper, the light
interpolator, the crown-depth quadrature — exists to *discover where to place nodes*
to hit a tolerance **without knowing the answer**. On a gradient pass the adaptive run
has already happened, so the node placement is already known. Re-discovering it is pure
overhead, and differentiating *through* the adaptive branching is fragile. So:

> **Record** where the double (adaptive) pass placed its nodes. **Replay** pinned to
> those nodes with the active scalar — no adaptive branching — and read the adjoints.

The double pass, once recorded, is **immutable**: its only job is to feed replays. This
is one idea. The apparent "levels" (§3) are the same idea at different depths, and they
collapse to a single System concept (§5).

---

## 2. The system, in four names

An AD run holds exactly three nouns and one verb; the whole design fits in them.

- **the double solver** (`d`) — the real adaptive run. R holds it. It is immutable
  after its pass, and it **owns the recording**.
- **the active solver** — the double System lifted to the active scalar (the *rebind*):
  the differentiable solver the gradient runs on. It holds **no** semantics between
  calls — it is re-seeded and re-fed the recording every call.
- **the recording** — the schedule plus the node stash the double pass produced, read
  *per call*.
- **replay** (the verb) — the active solver re-running the double's recorded schedule
  with the active scalar.

Reuse of the active solver and its tape (§7) is a *property* of it, not a concept of
its own: it is kept on the double solver and reused so an optimiser loop does not
rebuild it. It never changes a number.

---

## 3. Ownership sets the layers

The **Solver** is the generic engine: it steps arbitrary systems, either adaptively
(discovering node positions) or on a fixed grid (replaying them). The **System** owns
whatever field it builds to compute its rates. That ownership split *is* the layer
ladder:

| Layer | Recorded thing | Cadence | Owner | Needs |
|---|---|---|---|---|
| **L1** | the step **schedule** (`recorded_steps()` = `times()`) | per accepted step | **Solver** | nothing — always present, always differentiable |
| **L2** | the System's adaptive **node positions** (interpolator knots, `QK`/QAG subdivision) | **per step** | **System** | *record positions* |
| **L3** | the System's **background field** values | **per RK stage** | **System** | *record positions and values, and freeze* |

L1 is universal: every ODE has a stepper, so every System replays its schedule via
`advance_fixed(recorded_steps())` with no System participation. L2 and L3 are the
System's own state; the stepper only signals cadence.

### 3.1 Positions are per step; values are per stage

The cadences follow from what each freeze is *for*.

**Positions freeze per step.** The only reason to freeze node positions is to keep the
parameter-dependent adaptive *branching* off the tape — the refinement decisions are
`double`-valued control flow. One representative position set per step removes the
branching completely, and it is representative *by construction*: the step-size
controller accepts a step only when the RKCK error estimate stays under tolerance —
only when the state, and hence the field's adaptive structure, varies within the step
by less than tolerance. Any residual difference between a step's representative
positions and a stage's ideal positions is below tolerance: the noise floor the
adaptive scheme already accepts. Recording positions per stage would be six times the
storage to resolve something already inside the tolerance band.

**Values freeze per stage.** A frozen background is different in kind: each RK stage
evaluates the field at a genuinely different state, and a frozen consumer must read
back the *exact* value that stage consumed — a real quantity, not a discretization
choice. A per-step value would be wrong, not merely coarse.

### 3.2 One recording, read as two slices

A single double pass records the **union** — positions *and* values — in one
accumulate/commit cycle (§5). Each consumer reads its slice:

| Consumer | Reads | On the active pass | Cross term |
|---|---|---|---|
| **resident / total** (live) | the recorded **positions** | recompute the field on frozen positions with active state → self-feedback flows | present |
| **mutant / invasion** (frozen) | the recorded **values** | read as `double` background — off the tape, no recompute; only the mutant's own state is active | zero |

The union is forced by cross-feedback: a `run_mutant` reads the *resident's* per-stage
values, so the resident's `run` must record values for mutants it cannot foresee. It
records them for free — the per-stage hook captures the positions (at stage 0) and the
values (each stage) together. plant records exactly this union (`stand_height_history`
positions and `environment_history` values, `patch.h`).

---

## 4. Generality: not the interpolator, not the environment

Because the recording is System state and the core hooks carry only cadence, both
freezes generalise past any one example:

- **L2 is any adaptive component.** The record hooks stash whatever *positions* a
  System's adaptive machinery chose — spline knots for an interpolator, the
  Gauss–Kronrod subdivision for a quadrature, or a *list* of position-sets for a System
  with several (a TF24 Patch records the light-spline knots and the crown-depth
  quadrature nodes through the identical hooks). Nothing in the concept or `derivs`
  names a knot.
- **L3 is any recomputable sub-state.** A frozen background is whatever doubles a System
  stashed per stage — a held-fixed sub-state, frozen state variables, a canopy — not
  specifically an interpolated field. `derivs` calls `set_ode_state(y, index)` and the
  System decides what that reads.

"Environment" is a plant term; the odelia-neutral word is **field** — the background
coupling a rate reads. odelia grows no `Recording` noun and no replayable-interpolator
class: the schedule is `times()` (Solver state), the positions and values are thin
System state, and the numeric (interpolator, quadrature) is stateless and rebuilt from
the recording.

---

## 5. One `Replayable` concept: three signals and one query

There is **one** System, driven live or frozen per call (`run` vs `run_mutant`), and
one concept the stepper dispatches on behind `if constexpr (Replayable<System>)`. An
absent hook makes every call site a zero-cost no-op; nothing forces a System to be
differentiable or replayable.

```cpp
template <class S>
concept Replayable = requires(S s, int stage) {
  s.record_stage(stage);                                   // accumulate (per RK stage)
  s.record_ode_step();                                     // commit    (per accepted step)
  s.replay_step();                                         // load      (per step, active pass)
  { s.has_recorded_field() } -> std::convertible_to<bool>; // the frozen-field query
};
```

The three signals are the two nested loops of the adaptive stepper, not two recordings:

| Signal | Fires | Job |
|---|---|---|
| `record_stage(k)` | per RK stage, record pass | accumulate the union (positions@0, value@k) into scratch |
| `record_ode_step()` | per **accepted** step, record pass | commit the scratch → recording |
| `replay_step()` | per step, active pass, before the stages | load this step's positions / advance the recording index |

`record_stage` and `record_ode_step` are two because adaptive step **rejection** forces
commit-on-accept: a rejected step re-runs its stages, harmlessly clobbering the scratch;
only an accepted step commits. Folding the commit into the per-stage hook would record
rejected steps. `replay_step` is separate because positions must load *before* a step's
stages run. Three signals, one query — the minimal seam. (The per-step commit takes the
`_ode_` infix because a System's `record_step()` already serializes one collected state
to a history row; the two must not collide.)

### 5.1 The runtime state is two bits

A Replayable System needs to know only two things at runtime: **is it recording?**, and
when replaying, **is the field frozen?** (mutant) or recomputed live (resident). Those
are the only distinctions the hooks and the query branch on:

- `record_stage` / `record_ode_step` act only while **recording**.
- `has_recorded_field()` reports whether the L3 field cache is populated — the query
  `derivs` reads to route mutant replay.
- `replay_step` acts whenever replaying, which is *derived*, not stored: a System is
  replaying exactly when it is **not** recording and **has** a recording to read.

The core sees only the four concept members; how a System backs them is its own business.

---

## 6. Control flow of an AD run

Three runs share the identical stepper, `advance_fixed`, and tape machinery. They differ
in exactly three slots: which **system** the active solver carries, which **slice** of
the recording it reads, and the `has_recorded_field()` **query** inside `derivs`.

### 6.1 `run` (double) — the record pass

```
d.advance_adaptive({0, T})                                   [recording]
  per adaptive step  step():
    retry loop:
      stepper.step():  6 RK stages
        per stage k:  derivs(sys, y, k, t):  set_ode_state(y, t)   ← field on ADAPTIVE positions
                      record_stage(k)                              → accumulate: positions@0, value@k
      accept? ─no─→ shrink h, undo y/t, retry   (scratch clobbered)
             └yes─→ record_ode_step()                              → COMMIT scratch → recording[step]
  ⇒ recording = schedule times() (L1) + positions/step (L2) + values/step×stage (L3)
     owned by d, immutable hereafter
```

### 6.2 `run` (active) — resident / live gradient

The active solver carries the resident system, reads its **own** recording, and
recomputes the field so self-feedback flows.

```
active = active_solver(d)                                    [replay, live field]
  schedule  → advance_fixed grid          (L1, Solver→Solver)
  recording → active.system  (positions read, values ignored) (L2, System→System)
  tape on; computeJacobian(trait/IC seeds, forward):
    forward(x):
      active.system.scatter(x, slots); active.reset()
      active.advance_fixed(recorded_steps()):   ← replay schedule, NO adaptivity
        per step step_to:
          replay_step()                    → load THIS step's frozen positions
          stepper.step(): per stage derivs(...,k):
            query false → set_ode_state(y, t)   REBUILD field on FROZEN positions,
                                                ACTIVE values (gradient flows; self-feedback)
      return state
    seed adjoints → sweep → ∂functional/∂trait   (with self-feedback)
  tape off
```

### 6.3 `run_mutant` (active) — mutant gradient

The active solver carries the mutant system, reads the **resident's** recording as fixed
background, and the field is read back off-tape (its derivative zero).

```
active = active_solver(d_resident)  (mutant system)          [replay, frozen field]
  schedule           → advance_fixed grid   (L1)
  RESIDENT recording → active.system  (values read, positions ignored)  (L3, System→System)
  tape on; computeJacobian(MUTANT seeds, forward):
    forward(x):
      active.system.scatter(x, slots); active.reset()  ← background ← recorded value[0]
      active.advance_fixed(recorded_steps()):
        per step step_to:
          replay_step()                    → advance index (positions unused)
          stepper.step(): per stage derivs(...,k):
            query true → set_ode_state(y, k)   read recorded value @stage k as DOUBLE
                                               background (off tape, ∂/∂field = 0)
      return mutant state
    seed adjoints → sweep → ∂functional/∂mutant-trait   (field contribution 0)
  tape off
```

The only forks are inside `derivs`, keyed by the query; the stepper never learns what a
position *is*. The differentiation trick, stated once: on replay the node **positions**
are frozen doubles (no adaptive branching to corrupt the tape) while the node **values**
stay active — so the gradient flows through *what the nodes hold*, never through *where
they sit*. L3 goes one step further and takes the field off the tape entirely.

---

## 7. Reuse the active solver, read the recording per call

A gradient call builds an active solver, records a tape, sweeps it, and reads adjoints.
An optimiser loop (or a batch of mutants against one recording) calls this repeatedly.
Rebuilding the active solver and reallocating the tape every call is waste; the reuse
stays invisible to R.

**The active solver is the only cached thing.** The gradient runs on it, and a `Solver`
carries its own `tape`, so that tape *is* the reused tape — there is nothing else to
cache. It is held on the double solver, and its type is **named, not erased**: the
System supplies `rebind`, so `System::rebind<active>` spells the active solver's type
from inside `Solver<System>` — no `void*`, no `static_cast`.

**Anchored on the Solver object, not an R handle.** plant holds the solver as a plain
C++ member (`scm.h`: `Solver<patch_type> solver;`) and never wraps it in an XPtr, so an
anchor on the R handle would be invisible to it. Anchoring the active solver on the
`Solver` object gives the SCM the reuse for free.

**The recording is read per call, never frozen into the active solver.** The schedule
and node stash live on the immutable double solver/System and are handed over on every
call — not snapshotted at first build, not carried through `rebind` (which lifts values
only), not smuggled onto the solver as fit state. Each consumer hands its own slice: the
calibration entry hands its observations (in the `least_squares` functional), a
record→replay System hands its recording. The number a gradient returns comes entirely
from the per-call recording and seeds.

**Validity domain.** The recording is keyed to the ICs and params of the double run;
those fix the schedule and the positions. A cached run may be replayed with a different
**mutant** or a different **functional / observations**. Changing ICs or params
invalidates the recording and forces a re-record — reading it per call is what makes
that pickup automatic rather than a stale reuse. `has_recording()` and
`recorded_steps()` are the read-only surface behind the "forgot to record" guard.

Two senses of "cache" stay distinct: **cache = the amortized active solver + tape
(speed); record/replay = the recording (semantics).** Reusing the cache never changes a
number; the recording is what does.

---

## 8. Calibration is one functional, not the foundation

The foundational blocks are **an arbitrary functional** — a **pure reduction**: it reads
a replayed System and returns the scalar(s); it does not drive the solver — and
**replay-fixed** (this doc): the driver replays the recorded schedule and hands the
functional a positioned solver. `least_squares` is one functional built on them: it owns
only its measured **observations** and which recorded steps they attach to, and reads the
model's state at those steps. The schedule is the recording's, never the functional's. An
emergent functional reads native state at the recorded steps and carries no observations.
Calibration is one case among many, not the primary one — and its data lives in the
functional, not the solver.

---

## 9. The odelia-native demonstrator

Lorenz and leaf_thermal have no adaptive sub-numerics, so `Replayable` is exercised by a
purpose-built System: the odelia-native shrink of FF16's resident light — a scalar state
whose rate reads a field built by an adaptive interpolator over state-dependent node
positions. It carries one differentiable input and exercises all three depths against
finite differences:

- **L1 — schedule freeze.** Record `times()` on the double solver; replay `advance_fixed`
  on the active solver. Value reproduces to floating point; gradient matches central FD.
- **L2 — interpolator freeze (live).** Record the knots per step; replay with frozen
  positions carrying active values, field recomputed. Value + gradient-vs-FD.
- **L3 — frozen field.** Read the recorded per-stage values as `double` background; the
  trajectory reproduces the resident and the field's gradient contribution is zero (the
  invasion property in miniature). The field is plain double data — never an active
  constant with a zeroed derivative.
- **Reuse.** On a persistent double solver, repeated gradient calls reuse the active
  solver and tape and reproduce; a recording from a *fresh* double pass is picked up per
  call (the anti-staleness property an L1-only cache cannot express); live and frozen
  replays share one active solver with the mode chosen per call; replay-before-record
  errors.

---

## 10. Mapping to plant (FF16 / TF24 / TF24f)

The Systems differ only in **which node-sets** they record — a *list* of position-sets,
each read live or frozen:

- **FF16** records the light-spline knots. Resident → live recompute; mutant → frozen.
  The direct analogue of §9.
- **TF24 / TF24f** additionally record the crown-depth quadrature nodes and route the
  leaf-optimizer sensitivity through a supplied derivative. TF24f's coupled resident
  feedback is the stiff long-horizon boundary.

The plant `QK<S>` quadrature is out of odelia scope, but it is the same shape as the
interpolator — record positions, replay fixed, choose frozen or live — so plant adds a
quadrature node-set with no new mechanism. The single-concept, positions-vs-values model
is what guarantees that.
