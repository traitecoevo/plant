# AD infrastructure roadmap (plant × odelia)

**What.** A design for computing exact **trait gradients** of the `plant` SCM's
emergent outputs — how stand properties (LAI, biomass, basal area, offspring
production) respond to plant traits — by **automatic differentiation (AD)**, built on
`odelia`'s AD-aware ODE runtime rather than a plant-private AD stack.

**Why.** A validated prototype exists
([traitecoevo/plant#553](https://github.com/traitecoevo/plant/pull/553), tracking
issue [#472](https://github.com/traitecoevo/plant/issues/472)) but rebuilt inside
plant machinery that odelia already provides. This roadmap reaches the same gradients
with **surgical changes** to the existing plant components and a **small generic
surface** added to odelia — smaller, more maintainable, lower technical debt.

**Status.** Design phase — no code changes yet. Branches
`claude/ad-infrastructure-design` are prepared on the plant and odelia forks for the
eventual implementation.

---

## Read in this order

1. **[`ad-infrastructure-design.md`](./ad-infrastructure-design.md)** — the design.
   The thesis, the settled decisions, the three layers (odelia / plant / UX), the
   surgical changes, the two workflows, and the four replay levels. Start here.
2. **[`ad-r-interface.md`](./ad-r-interface.md)** — the R/C++ boundary: why XAD types
   are awkward through Rcpp, the "only doubles cross" invariant, and **user stories**
   for each persona. Read after the design.
3. **[`ad-issues.md`](./ad-issues.md)** — the work breakdown: scoped items
   (odelia / plant / R-boundary / prototypes), dependencies, and a critical-path
   build order. Read when planning implementation.

---

## Key concepts (glossary)

For a reader new to the model or to AD. Fuller treatment is in the design doc.

- **SCM** — the plant "Solver for Characteristics Method": integrates a size- and
  patch-structured population as cohorts introduced on a schedule and stepped by an
  ODE solver. The **Patch** is the state being integrated; it is already an
  `odelia::ode::Solver` System.
- **Emergent output / metric** — a stand-level property that emerges from the cohorts:
  LAI, biomass, basal area, offspring production. These are the quantities we
  differentiate.
- **Trait gradient** — the derivative of an emergent metric with respect to plant
  traits (lma, wood density, …). The deliverable.
- **Two workflows.**
  - *Resident / total* — the whole stand differentiates **with self-feedback**: a
    trait change re-shades the canopy the stand grows in.
  - *Mutant / invasion* — a **rare** mutant's fitness gradient against an established
    resident whose canopy is held fixed (the mutant is too rare to shade itself).
    This is the selection gradient of evolutionary ecology.
- **Birth rate & demographic equilibrium** — beyond trait gradients, the derivative of
  the net reproduction ratio with respect to birth rate (`d R0/d birth_rate`) drives the
  Newton solve for the equilibrium birth rate (where R0 = 1) — the density at which a
  strategy sustains itself. The birth-rate gradient is a first-class deliverable, not
  only a sensitivity.
- **Reverse-mode AD / XAD / tape** — reverse-mode records operations on a *tape*, then
  sweeps it backward to get all input derivatives from one output. **XAD** is the AD
  library odelia vendors and compiles once. Reverse mode is optimal here because there
  are many traits (inputs) and few metrics (outputs).
- **Replay levels (L0–L3)** — the SCM has several *adaptive* constructions (the node
  schedule, the ODE step sizes, the light interpolator, the crown quadrature). AD
  needs each frozen to its recorded placement so the result is differentiable. The
  key idea (design §7.5): run once adaptively, **record where the nodes landed, replay
  on them fixed** with the active scalar.
- **odelia** — the family's AD-aware ODE runtime: compiles XAD once, ships a
  scalar-templated `Solver`, a reverse-mode gradient driver, and a differentiable
  spline. plant already `LinkingTo: odelia`.
- **The spike** — PR #553, the validated prototype. In this roadmap it is the
  **specification and the regression oracle**, not the code that ships.

---

## The one-sentence design

Make the existing plant Patch and Strategy AD-compatible with small in-place changes,
differentiate the SCM runs plant already has (`run` for resident, `run_mutant` for
invasion) on a recorded-then-fixed replay, evaluate a plant-supplied emergent
functional, and let odelia own the generic AD mechanism (tape, Jacobian, functional
seam) it already almost provides — so only `double` ever crosses back to R.
