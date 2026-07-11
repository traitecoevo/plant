# Multi-patch water redistribution: architecture proposal (#427)

Status: proposal / design. Branch `feature/multipatch-water-proposal`, based on
`develop`. No implementation yet — this is the design doc for the water-flow arm
of the multi-patch epic
[#427](https://github.com/traitecoevo/plant/issues/427), to be linked from that
epic's "## Issues" list (currently "to be linked").

> **Scope.** This realises #427's "patches connected by water flow between them"
> acceptance criterion and its graph formulation (nodes = patches, edges =
> flows). It builds on the TF24 hydraulics line
> ([#424](https://github.com/traitecoevo/plant/issues/424)) and the per-patch
> N-layer soil-water store now on `develop`, and it depends on the discrete-events
> mechanism scoped in
> [#522](https://github.com/traitecoevo/plant/issues/522). See *Relationship to
> existing issues* at the end for the full wiring.

## Motivation

`develop` now carries a per-patch N-layer vertical soil-water store (currently
N = 15) with a Sperry-style hydraulic link to carbon gain (FF16w line). The
natural next step is letting patches exchange water so we can model spatially
structured water limitation — the semi-arid, event-pulsed regime (Mulga
banding, run-on/run-off) being the primary target, with a design that
generalises in principle across the aridity gradient.

This issue proposes an architecture that adds spatial coupling **without**
introducing a stiff coupled ODE system across patches, staying consistent with
our "design out stiffness rather than reach for implicit solvers" philosophy
(cf. [#529](https://github.com/traitecoevo/plant/issues/529), taming TF24f
fast-relaxation stiffness locally rather than with a global implicit solver).

## Core idea: three separable water processes

Rather than one lumped hydrology, the design treats water movement as three
physically distinct, independently-toggleable processes. Which are active is
set by climate/soil, not by a hand-placed regime switch:

1. **Surface lateral redistribution at events** — between-patch, algebraic
   reallocation at each discrete rainfall event. Dominant in drylands
   (infiltration-excess run-on/run-off). Weights are state-dependent
   (cover/infiltration contrast) + topographic.
2. **Vertical infiltration + drawdown within a patch** — the existing N-layer
   bucket (N = 15). Always active. This is where the Sperry link and the
   dendro-relevant multi-year deep-store buffer live. *Within-patch only — not
   a coupling term.*
3. **Lateral subsurface redistribution between patches** — continuous flow via
   a water table. Effectively off in semi-arid systems (emerges from K(θ): low
   soil moisture → collapsed unsaturated conductivity → negligible lateral
   flow, no hand toggle needed), on in wetter systems.

Mulga: 1 + 2 on, 3 off. Wet forest: 2 + 3 on, 1 off. Sub-humid: all three.

## Proposed implementation

### 1. A `Landscape` / `MetaPatch` container

Introduce a container above `Patch` holding `N` patches plus topographic
metadata. This is where cross-patch coupling lives; `Patch` and `Environment`
stay largely as-is so the demographic/SCM machinery is untouched. How this
container is *presented to users* (constructing a landscape, running it) should
be scoped alongside the unified run-interface work in
[#428](https://github.com/traitecoevo/plant/issues/428).

Per-patch static attributes (precomputed once):
- wetness/run-on index (`ln(a/tanβ)`-style) → sets redistribution weights
- **heat-load index** (f(slope, aspect, latitude)) → scalar multiplier on PET /
  evaporative demand. Carries the aspect story (hot equator-facing ridge vs
  cool pole-facing gully) through a single per-patch scalar that then cascades
  into θ and hence into process 3 via K(θ).
- **soil profile properties** — see §5, these covary with topographic position.
- topographic order (for the downslope sweep, below)

### 2. Event-time redistribution operator (process 1)

Driven by a discrete rainfall event schedule. This depends on the discrete-events
mechanism proposed in [#522](https://github.com/traitecoevo/plant/issues/522)
(event-driven hydraulics: integrate to the event, apply an algebraic jump to soil
water, resume) — which itself generalises the existing node-introduction
machinery into a first-class `(time, action)` events queue. Between-patch
redistribution is a new event *action* on that queue.

At each event:
- distribute incident water across patches via a **state-dependent**
  redistribution operator: weights combine topographic term + cover/
  infiltration term. On flat ground the state-dependent term dominates →
  banding can self-organise; on steep ground the topographic term dominates →
  externally imposed flow direction. Same operator, regime set by topography.
- implemented as a **single downslope sweep** in topographic order (cascade,
  can re-infiltrate downslope). Algebraic, O(N), differentiable.
- closed budget for small/medium events; spill/deep-drainage loss term for
  large events.

Keep the operator **smooth** (no hard infiltration thresholds) to preserve
clean AD gradients — important for downstream GEK emulation / calibration
(the AD programme is [#472](https://github.com/traitecoevo/plant/issues/472),
with the finite-difference inventory in
[#537](https://github.com/traitecoevo/plant/issues/537); any new derivative
introduced here should be AD-able rather than FD).

### 3. Between-event lateral flux (process 3): donor-controlled

For the wet-region generalisation, lateral subsurface flux is made
**donor-controlled**: flux out of a patch depends only on the *upslope* patch's
storage, not the head *difference*. This makes the coupling a DAG (strictly
downslope), so the Jacobian is **triangular** in topographic order → solvable by
a single downslope sweep, no coupled implicit solve, AD-friendly.

Note this triangularises coupling but doesn't remove fast *local* drainage: near
saturation K(θ) is large → stiff diagonal. But that's now a *decoupled* per-patch
fast timescale, which our QSSA/relaxation toolkit
([#529](https://github.com/traitecoevo/plant/issues/529)) handles directly. K(θ)
self-diagnoses which patches are stiff (wet ones) vs not (dry Mulga patches →
not stiff, semi-arid case stays trivial).

Flux magnitude via a relaxation timescale τ:
- τ → ∞ : decoupled, event-redistribution-only (semi-arid; process 3 off)
- τ → 0 : instant equilibration to topographic index (wet-region / TOPMODEL limit)
- intermediate : sub-humid

Donor-control caveat: accurate in steep/gravity-dominated terrain; first-order on
flat/near-saturated ground where the receiver's state genuinely backs up flow.
Acceptable for "generalises in principle"; flag as the loosest approximation.

### 4. Resource hand-off (why this generalises ecologically)

The water machinery only does ecological work where water is limiting. Toward
the wet end, water stops being limiting and the existing light-competition
kernel should dominate. The general model is **water routing + light competition
running simultaneously**, aridity determining which binds. This is a
co-limitation composition, not a τ knob — worth designing the Environment so
both kernels compose cleanly. The cross-patch water field is the water analogue
of the light kernel, but with an upslope→downslope *asymmetry* light competition
lacks — flagging because it may interact with symmetry assumptions in the
fitness machinery.

### 5. Topographically covarying soil properties

Soil depth, texture, water-holding capacity and the K(θ) curve are **not**
independent of topographic position — they covary with it systematically along
the catena:
- downslope/convergent: colluvial, deeper, finer, higher capacity
- upslope/divergent ridge: erosional, shallower, stonier, lower capacity

This **reinforces** the redistribution gradient rather than being orthogonal to
it: the downslope patch receives more water *and* has more capacity to hold it
*and* (being wetter, higher θ) has higher K(θ) so conducts it onward. The
gully-as-refugium story is a triple stack — cool aspect + run-on position + deep
soil — not just the first two.

Implementation:
- soil properties become **per-patch static attributes** (layer capacities,
  retention/K(θ) parameters), not a shared profile with only forcing differing.
- recommend **fixed N (=15) across patches with varying layer thickness / depth**,
  so the per-patch ODE state vector keeps a uniform shape (simpler for the
  stepper) while a shallow ridge profile is represented by thinner layers than a
  deep gully profile.
- benign interaction with the stiffness diagnosis: downslope patches are both
  wetter and finer → the stiff (high-K, fast-draining) ones → and they sit
  *last* in the downslope sweep, after their inputs are determined. Stiff
  patches are encountered at the end of the topological order, not the start.

**Inference hazard (flag now):** because soil capacity covaries with the
redistribution weights, a wet downslope patch can be explained either by "receives
run-on" or by "high water-holding capacity" — the redistribution parameter and the
soil-capacity parameter may be only *jointly* identifiable, not separately.
Mitigate with informative priors on soil properties from independent soil data /
partial pooling (hierarchical) rather than trying to estimate both freely from
moisture or growth data alone. Much cheaper to design for than to discover during
calibration.

## Stepper structure: mirror the existing outer/inner loops per-patch

> Terminology note: `NodeSchedule` is being renamed to `Schedule`
> ([#505](https://github.com/traitecoevo/plant/issues/505)); this section uses
> the current `build_schedule` / node-introduction names.

The current SCM uses **two adaptivities at different levels**, not a single
unified stepper:

- **Inner, online — adaptive time stepping.** During a single forward
  integration the adaptive ODE stepper (`ode_tol_rel`/`ode_tol_abs`,
  `ode_step_size_max`) reacts to local stiffness/curvature, with node
  introduction times as fixed breakpoints it must land on. The node schedule is
  *fixed* for the duration of a run.
- **Outer, offline — adaptive node refinement.** `build_schedule` runs the
  entire patch to completion, assesses **cohort-spacing quadrature error** (a
  size-density integration error via `run_scm_error` / `area_leaf_above`, *not*
  the ODE error), inserts introduction times where error > `schedule_eps`, and
  re-runs the whole trajectory. Repeats to consistency.

Two different error notions, two different loops. The multi-patch design should
**mirror this two-level structure at the landscape level**, and both loops
sweep in topological (downslope) order courtesy of donor-control:

- **Outer (node refinement) across coupled patches.** Because patch B's
  environment depends on patch A's water-uptake trajectory, schedules **cannot**
  be refined independently per patch. But donor-control (one-way downslope
  coupling, §3) makes the refinement *sweepable in topological order*: fully
  refine the most-upslope patch (no water input from others), freeze it, refine
  the next patch down against that frozen upstream trajectory, and so on. The
  same DAG that makes the flux Jacobian a downslope sweep makes schedule
  refinement a downslope sweep — a real payoff of donor-control beyond the flux
  solve. Without donor-control this becomes a joint refinement over all patches'
  schedules (an outer-outer landscape-consistency loop) — much more expensive.
- **Inner (time stepping) — global vs per-patch, set by whether process 3 is
  on.**
  - Semi-arid (process 3 off): between events patches are genuinely decoupled,
    so run **per-patch independent adaptive steppers synchronised only at
    rainfall-event breakpoints**. A dry, slow patch is not dragged onto the tiny
    time steps a wet, fast neighbour demands. Correct *and* more efficient.
  - Wet / sub-humid (process 3 on): continuous coupling wants a common time
    axis → either a **global stepper** over the concatenated all-patch ODE
    state, or a multi-rate scheme (more work). Start with the global stepper.

This is the part most likely to bite whoever implements it, so it's called out
explicitly rather than left to discovery.

## Suggested sequencing

1. `Landscape`/`MetaPatch` container + per-patch static attributes (heat-load,
   wetness index, **per-patch soil profile**), no coupling yet — patches
   independent. Validates the container doesn't disturb SCM/fitness, and that
   per-patch heterogeneous soil works with the existing N-layer bucket.
2. Event-time surface redistribution (process 1), closed budget, per-patch
   steppers synced at event breakpoints. Reproduce prescribed-flow (steep) case
   first, then emergent banding (flat). *(Requires the #522 events queue.)*
3. Donor-controlled between-event flux (process 3) + per-patch QSSA (#529) +
   global inner stepper. Unlocks sub-humid and wet regimes.
4. Simultaneous light+water limitation composition.

Steps 1–2 alone cover the semi-arid target and are cheap (algebraic coupling
only, per-patch independent steppers). 3–4 are the generalisation and can wait.

## Validation checks

- Flat-ground emergent banding: reconcile emergent wavelength + band migration
  rate against Klausmeier/Rietkerk predictions. Note our nonlocality is
  *derived* from explicit water routing rather than a phenomenological
  competition kernel — a genuine point of novelty worth being explicit about.
- Steep terrain: spatial trait sorting along the fixed upslope→downslope fitness
  gradient.
- Semi-arid limit recovers current single-patch FF16w behaviour when N = 1.

## Open questions

- **Extend the existing SCM or build a separate multi-patch SCM?** The core
  build decision. Two ends of the spectrum:
  - *Extend in place* — teach the current `run_scm` / `Patch` / `SCM` machinery
    to hold `N` patches (single-patch = the N = 1 special case). Pro: one code
    path, so the semi-arid N = 1 limit recovers current FF16w behaviour *by
    construction* rather than as a validation target, and no duplicated
    demography/fitness logic to keep in sync. Con: risks entangling the
    landscape-coupling concerns into hot single-patch paths (the `[speed]` work,
    e.g. #471) and the fitness/equilibrium machinery that currently assumes one
    patch.
  - *Separate `Landscape`/`MetaPatch` layer* — a new container that *composes*
    unmodified `Patch`/`SCM` objects and owns only the cross-patch coupling
    (redistribution operator, downslope sweep, landscape-level schedule
    refinement). Pro: keeps `Patch`/`Environment`/fitness untouched (the stated
    §1 goal), isolates the new stiffness/stepper concerns, and matches the
    graph-of-nodes framing in #427. Con: the two-level stepper (per-patch vs
    global inner loop, §"Stepper structure") straddles the boundary, so the
    container can't be a pure black-box wrapper — it needs to reach into
    per-patch stepping when process 3 is on.
  - Leaning: **composition (separate layer) as the default**, because the
    semi-arid target (steps 1–2) is purely algebraic coupling over independent
    per-patch steppers — exactly what a thin container gives cheaply — and it
    keeps the single-patch fast path clean. The extend-in-place pressure only
    arrives with process 3 (global inner stepper over concatenated state), which
    is step 3. Decide before step 1, since it sets where the state vector and the
    stepper ownership live. Cross-cuts the unified-interface scoping in #428.
- Does `Environment` need restructuring to carry two simultaneous limiting
  resources cleanly, or can water enter purely via the strategy's PET/Sperry
  path with light unchanged?
- Where should the rainfall-event schedule live — extend `Disturbance`/`Control`,
  or a new driver object? This should be resolved *with*
  [#522](https://github.com/traitecoevo/plant/issues/522) (the events queue is
  proposed to be "passed in alongside or as part of the drivers") and
  [#425](https://github.com/traitecoevo/plant/issues/425) (drive model with
  variable environments); note the driver-iteration refactor #364 (closed)
  already moved env drivers to a single list-handling chunk, which is the natural
  place to hang an events list.
- Fitness machinery + upslope→downslope asymmetry: does the directed spatial
  structure break any symmetry assumptions in the current invasion-fitness /
  equilibrium solver?
- **Building the patch-graph for an arbitrary landscape.** The whole scheme
  assumes a topologically-ordered DAG of patches (donor-control, downslope
  sweep). We need a general way to *construct* that graph from a landscape
  description. This is the concrete realisation of #427's graph formulation
  (nodes = patches, edges = flows); the spatial-grid variant #427 also asks for
  is the lattice special case. Open sub-questions:
  - Input: a DEM (raster), an idealised synthetic hillslope, or an abstract
    adjacency the user supplies? Support which, first?
  - Flow-direction / accumulation algorithm to turn a DEM into edges: single-flow
    (D8 — gives a clean tree, one downslope neighbour per patch, trivially
    sweepable) vs multiple-flow (D∞/MFD — more realistic divergent flow on
    ridges/fans, but a patch then has several downslope receivers with split
    weights; still a DAG, still sweepable, but edges carry fractions). D8 is the
    simpler starting point; MFD is probably needed for realistic fans/banding.
  - **Flat ground has no unique flow direction** — exactly the emergent-banding
    regime. On truly flat terrain the DAG is degenerate/ill-defined, so the graph
    can't come from topography alone; edges must be induced by the
    state-dependent (cover/infiltration) term and can *reorganise* as vegetation
    changes. Implies the graph may need to be **dynamic** (recomputed as cover
    evolves) in the flat/emergent case, vs **static** (topographically fixed) in
    the steep case. Decide whether to support both from the outset or ship
    static-first.
  - Pit/depression handling and the guarantee of acyclicity: standard DEM
    pre-processing (fill/breach sinks) ensures no cycles → keeps the sweep valid.
    Need to confirm whatever we adopt cannot produce a cycle, or the triangular-
    Jacobian argument fails.
  - Patch = raster cell, or patch = an aggregated hydrological unit (hillslope
    element / catena position)? Aggregation keeps N tractable and matches the
    demographic scale better than raw cells.
- **How do we verify this is working reasonably?** Beyond the validation checks
  above, a graded ladder of tests:
  - *Conservation / sanity:* closed water budget across the landscape between
    events (in = stored + transpired + spilled), to machine tolerance for the
    algebraic reallocation. First and cheapest check. (This is the landscape-scale
    analogue of the single-patch conservation check already tracked in
    [#454](https://github.com/traitecoevo/plant/issues/454).)
  - *Analytic limits:* N = 1 recovers single-patch FF16w exactly; τ → ∞ recovers
    fully-decoupled patches; τ → 0 recovers the steady-state topographic-index
    (TOPMODEL) water-table shape. Each is a known target.
  - *Idealised hillslope:* a monotonic slope should produce a monotonic
    upslope→downslope moisture/biomass gradient with no spurious structure;
    check moisture profile against a hand-solvable donor-control cascade.
  - *Emergent pattern:* on flat ground, does banding wavelength + migration rate
    fall in the range predicted by Klausmeier/Rietkerk for comparable
    parameters? (Not exact — our feedback is mechanistic, not phenomenological —
    but should be in the right ballpark and scale correctly with rainfall.)
  - *Convergence:* results stable under patch-graph refinement (halving cell
    size / doubling N shouldn't move emergent properties much) and under
    `schedule_eps` / ODE-tolerance tightening — i.e. numerical vs ecological
    signal are separable.
  - *Empirical, if reachable:* Mulga grove/intergrove moisture and growth
    contrasts, or the dendro series, as a weak external check — noting the
    soil-vs-redistribution identifiability caveat (§5) means this validates the
    *joint* system, not individual parameters.
- Node-schedule refinement across coupled patches: confirm the topological-order
  sweep converges in practice (does freezing an upstream schedule while refining
  downstream ever need a second pass once downstream refinement shifts nothing
  upstream? It shouldn't, given one-way coupling — but worth a test).
- Fixed-N / varying-thickness soil layers: confirm the Sperry hydraulics and
  K(θ) behave sensibly when layer thickness differs markedly between a shallow
  ridge and a deep gully patch.

## Relationship to existing issues

**Parent**
- [#427](https://github.com/traitecoevo/plant/issues/427) **[patch variations]
  Multi-patch dynamics** — this proposal is the water-flow arm of that epic; the
  graph formulation and the "connected by water flow" acceptance criterion are
  realised here. Link this doc under #427's "## Issues" (currently "to be linked").

**Depends on**
- [#522](https://github.com/traitecoevo/plant/issues/522) **Enable event-driven
  hydraulics** — the discrete `(time, action)` events queue that process 1's
  event-time redistribution runs on. Hard dependency for step 2.
- [#424](https://github.com/traitecoevo/plant/issues/424) **new growth model with
  hydraulics and competition for water** — the TF24 line that introduced the
  per-patch N-layer soil-water store + Sperry link this builds on.
- [#529](https://github.com/traitecoevo/plant/issues/529) **Tame TF24f
  fast-relaxation stiffness (IMEX / operator-split / QSSA)** — the per-patch
  fast-timescale toolkit process 3 delegates decoupled stiffness to.

**Interacts with**
- [#425](https://github.com/traitecoevo/plant/issues/425) **Drive model with
  variable environments** and #364 (closed) **driver iteration** — where the
  rainfall-event schedule and per-patch drivers plug in.
- [#428](https://github.com/traitecoevo/plant/issues/428) **scope user interface
  for running… multiple patches** — the user-facing surface for
  `Landscape`/`MetaPatch`.
- [#472](https://github.com/traitecoevo/plant/issues/472) **[AutoDiff] epic** /
  [#537](https://github.com/traitecoevo/plant/issues/537) **FD-derivative
  inventory** — the "smooth operator / clean AD gradients / GEK emulation"
  requirement; new derivatives here should be AD-able.
- [#454](https://github.com/traitecoevo/plant/issues/454) **Confirm water is
  conserved** — single-patch conservation check whose landscape-scale analogue is
  the first validation test above.
- [#505](https://github.com/traitecoevo/plant/issues/505) **Rename NodeSchedule →
  Schedule** — naming in the stepper section.

**Sibling `[patch variations]` run-modes** (not dependencies, but the same
container should accommodate them per #428):
[#519](https://github.com/traitecoevo/plant/issues/519) continuous patch,
[#520](https://github.com/traitecoevo/plant/issues/520) LSM mode.
