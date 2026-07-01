# Seed: TF24 NSC storage to buffer short-term productivity variation (#517)

Handoff brief for starting the NSC work in a fresh session. Everything below was
established while investigating the SCM cohort-density blow-up (#550). Read this
first so you don't re-derive it.

## Goal (#517, epic)

Add a non-structural-carbohydrate (NSC) **storage pool** to the TF24 model so that
growth and mortality respond to *buffered* carbon availability rather than
*instantaneous* net production. Motivation: currently a productivity drop
translates instantly into a mortality spike. Fiona Robinson's thesis showed this
is unrealistic — real plants draw down stored carbon, so they die gradually, not
in a sudden rush.

## Why this is the root-cause fix for #550

The SCM blow-up (#550) under extreme seasonal drought was traced to **mortality
stiffness**, not (as first hypothesised in the now-closed #551) a growth-gradient
discontinuity. Confirmed by direct instrumentation:

- `mortality_growth_dependent_dt = a_dG1 * exp(-a_dG2 * productivity_area)`
  reaches ~1e32+ under deep carbon deficit (`src/tf24_strategy.cpp:641`). At the
  blow-up: `dg/dh ≈ 0`, `mortality ≈ 8e32`.
- That astronomically stiff rate overflows the coupled `(log_density, mortality)`
  ODE under *any* explicit integrator (adaptive RKCK, tight tol, fixed Euler all
  fail identically).
- The growth-rate gradient `dg/dh` is smooth (~0.22) in height — no
  discontinuity. (#551 hypothesis refuted.)
- Related sharpness: growth has a hard cutoff `net_mass_production_dt > 0 ? grow
  : 0` (`src/tf24_strategy.cpp:186`), and the hydraulic optimum `opt_psi_stem`
  clamps to `psi_crit` as soil dries through θ≈0.16 (a KKT corner). Both are
  crossed sharply when carbon balance crashes.

### Local fixes tried and rejected (do NOT re-attempt as the fix)

All three shift *which* (mpl, amp, param) combinations blow up without robustly
preventing it — chaotic knife-edge. Evidence + tables are in #550 comments.

- Growth-cutoff smoothing (logistic/softplus on net production): moved the
  knife-edge, didn't remove it.
- Mortality cap (`min(..., a_dG3)`): chaotic in both cap value (mpl=30: 100,150,500
  pass; 50,200,300 fail) and drought amplitude (cap=100: amp 0.34–0.38 fail, 0.30/
  0.40 pass). Masks, doesn't resolve. Also a hard `min` is an AD kink (see below).
- Density-rate limiter (cap `d(log n)/dt`): blow-ups persist; soil-coupling needs
  `log_density ≲ 7`, which would flatten the distribution.

Conclusion: no local parameter tweak works. The blow-up is a coupled,
chaotically-sensitive density–soil–hydraulic instability near carbon shutdown.
NSC buffering (this issue) is the principled fix.

### Interim already shipped

`#552` (merged to `develop`, squash `cc4c76e1`): a graceful-failure guard
(`Patch::check_finite_ode_state`, `inst/include/plant/patch.h`) that stops with a
descriptive message when a cohort density or environment state goes non-finite.
It does NOT make blow-up regimes complete — it just fails clearly. Keep it; NSC
should make those regimes actually complete, at which point the #550 regression
test (`tests/testthat/test-strategy-tf24.R`, "SCM cohort-density blow-up fails
gracefully (#550)") must be updated (it currently asserts the error).

## Design direction for NSC (proposal — confirm with Daniel)

- Add an NSC **storage** state variable to TF24 (one extra ODE state).
- Net production charges/discharges storage; **growth and mortality read a
  smoothed "available carbon" signal from storage**, not instantaneous
  `net_mass_production_dt`. This buffers troughs → mortality no longer spikes to
  1e32, growth cutoff no longer hit sharply.
- Allocation priority to think through: maintenance/respiration → storage
  refill → growth → reproduction (and what happens as storage depletes).
- Mortality: tie growth-dependent mortality to **storage depletion / low reserves**
  rather than instantaneous productivity, so death is gradual as reserves run out.
- Parameters (new `TF24_Pars` fields): storage capacity (likely scales with leaf
  area / live mass), charge/discharge or turnover rates, and the
  reserve→mortality relationship. Calibration cues in Fiona Robinson's thesis.

### Prefer SMOOTH formulations (AD-readiness)

plant is moving to reverse-mode AD via XAD, through templated kernels
(`inst/include/plant/models/ff16_production_kernel.h`; #472/#537, currently FF16
net-mass-production w.r.t. traits; TF24 and mortality not yet in a kernel).
Fitness gradients will eventually pull mortality into AD. So:

- Use smooth saturations (softplus / logistic / logsumexp), NOT `min`/`max`/hard
  `?:` switches, in any carbon→growth/mortality mapping. A hard `min` compiles
  under XAD only if the constant is the templated scalar type, and it introduces
  a kink (discontinuous derivative, zero gradient in the capped region).
- The existing `net > 0 ? grow : 0` growth cutoff is itself an AD hazard worth
  smoothing as part of this work.

## Implementation touchpoints

- `src/tf24_strategy.cpp`: `compute_rates` (add storage rate; route growth/
  fecundity/heartwood and mortality through buffered carbon), `mortality_*`,
  `net_mass_production_dt`.
- `inst/include/plant/models/tf24_strategy.h`: `TF24_Pars` (new params);
  `state_size()` (+1 for storage); `state_names()`; index caching in
  `refresh_indices()`.
- `inst/RcppR6_classes.yml`: add new `TF24_Pars` fields + storage state → then
  regenerate: `make RcppR6` (RcppR6 is installed) and rebuild.
- Tests: `tests/testthat/test-strategy-tf24.R` (Defaults enumerates pars; Critical
  Names enumerates state names — both need updating; update the #550 graceful-
  failure regression test since those regimes should now complete; E-conservation).

## Validation targets

- The #550 repro sweep must **complete with bounded, sane densities AND be robust**
  across `mpl ∈ 20–70` and `amp ∈ 0.30–0.40` (unlike the mortality cap, which was
  chaotic across exactly this sweep).
- Healthy runs (`mpl=10/13`) unchanged or changed defensibly.
- Mortality rate stays finite/O(1–10) through drought troughs (no 1e32 spike).

## Reusable instrumentation harness

Isolated single-individual probe (cracked the earlier "insufficient number of
points" blocker — it was the uninitialised light environment):

```r
e <- Environment("TF24")
e$set_soil_number_of_depths(5)
e$set_soil_water_state(rep(theta, 5))
e$set_fixed_environment(1.0, 40)      # full light up to 40 m
pl <- TF24_Individual()
pl$set_state("height", h)
pl$compute_rates(e)
pl$aux("net_mass_production_dt"); pl$aux("opt_psi_stem"); pl$rate("height")
```

`#550` repro: `scm_base_parameters("TF24","TF24_Env")`, `max_patch_lifetime=mpl`,
`add_strategies(., trait_matrix(0.07,"lma"))`, 5 soil depths at 0.2, rainfall
`amp*sin(2*pi*x)+0.5` on `x=seq(0,mpl,length.out=mpl*6)`, then `run_scm`.

Build: `pkgload::load_all(".")` (~1–2 min after header/cpp change). TF24 full SCM
runs are slow (~min) but the blow-up fails fast (t≈3.75, first deep trough).

## Refs
- #517 (this epic) · #550 (blow-up bug, mechanism) · #551 (closed, refuted) ·
  #552 (merged graceful guard) · #472/#537 (XAD AD kernel roadmap).
