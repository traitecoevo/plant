# Penman-Monteith leaf energy balance in TF24 / TF24f

## Context

`penman_monteith_implementation.md` proposes replacing the simple Fick's-law
transpiration closure (`E = gs·VPD`, leaf temperature assumed = air temperature)
with a Penman-Monteith energy balance that solves latent heat flux and leaf
temperature simultaneously, then feeds the computed `Tleaf` into the Farquhar
photosynthesis temperature scaling. The motivation is that under hot, high-
radiation Australian conditions `Tleaf` can deviate 5-10 °C from `Tair`, biasing
both the water cost `E` and the carbon gain `A`.

This plan records **where PM would wire into TF24/TF24f**, and answers the two
questions the investigation asked: (1) does it change the complexity of the leaf
photosynthesis modules, and (2) what additional environmental variables are
needed. The user asked for an investigation for implementation — this is an
implementation plan plus a scoped complexity assessment, not a committed build.

**Key correction to the doc's framing.** The doc presents PM as *predicting* `E`
from stomatal resistance (`λE = (Δ·Rn + ρcp·VPD/ra)/(Δ + γ(1+rs/ra))`). In TF24
that would be wrong: `E` is **not free** — at every candidate operating point it
is already pinned by the Sperry hydraulic supply
(`E = transpiration(psi_stem, psi)`, `leaf_model.cpp:1126`), a purely hydraulic
quantity. Using the full PM inversion would *replace* that supply and gut the
demand-supply structure the model is built on. The correct integration uses only
the doc's **second** equation — the energy-balance residual
`Tleaf = Tair + (Rn − λE)·ra/ρcp` — with the hydraulically-known `E`. `Tleaf` is
then an **explicit algebraic function of `E`**: no PM inversion, no fixed-point
iteration, and (crucially) `A(Tleaf)` does **not** feed back into `E`, so there
is no `Tleaf↔A↔gs` loop within an operating point. This is what makes the feature
tractable.

The doc itself (Section 7) is explicit that PM should be a **staged hypothesis
test**: each component (PM core, Rn submodel, wind/`ra` profile) is kept only if
a sensitivity analysis shows it materially changes outputs. That framing should
be honoured — build the minimal core first, measure, then decide on the rest.

## How TF24/TF24f actually work today (relevant facts)

- TF24 and TF24f share one embedded `Leaf` submodel
  ([leaf_model.h](../../GitHub/plant-family/plant-dev2/inst/include/plant/leaf_model.h),
  [leaf_model.cpp](../../GitHub/plant-family/plant-dev2/src/leaf_model.cpp)):
  coupled Farquhar (FvCB) photosynthesis + Sperry-style hydraulic
  profit-maximisation. TF24f (`src/tf24f_strategy.cpp`) subclasses TF24 and
  replaces the per-step re-optimisation with a tracked ODE state advanced by
  gradient ascent on profit.
- **TF24 does not use Fick's law as its top-level closure.** Transpiration is the
  *hydraulic supply* through soil→root→stem vulnerability curves; the operating
  point is found by optimising the root-collar water potential
  (`find_root_collar_psi`, `leaf_model.cpp:798`). The Fick relation the PM doc
  targets is present as the **E↔gs conversion** in `Leaf::stom_cond_CO2`
  (`leaf_model.cpp:1164-1166`):
  `gs = atm_kpa·E·kg_to_mol_h2o / atm_vpd / 1.67`.
- **Leaf temperature and VPD are prescribed constant extrinsic drivers**, not
  solved. Registered in the `TF24_Environment` constructor
  (`tf24_environment.h:82-88`: `leaf_temp`=25, `atm_vpd`=1, plus `PPFD`,
  `rainfall`, `ca`, `atm_o2_kpa`, `atm_kpa`), read via `get_leaf_temp()` /
  `get_atm_vpd()` and passed into `Leaf::set_physiology(...)`
  (`leaf_model.cpp:175`) at the single call site `tf24_strategy.cpp:440`.
- The Farquhar biochemistry is **already fully temperature-dependent**:
  `vcmax_`, `jmax_`, `gamma_`, `kc_`, `ko_`, `R_d_` are Arrhenius/peaked-Arrhenius
  functions of `leaf_temp_` (`leaf_model.cpp:246-253`). Because `leaf_temp_` is
  constant, that block is computed **once per `set_physiology`** and memoised
  (`photo_temp_cached_`).
- Absorbed PAR is already available per crown light point: `radiation = k_I ·
  light · PPFD` (`tf24_strategy.cpp:447-450`), so shortwave `Rn` needs no new
  driver.
- TF24f's default gradient is an **exact analytic AD/implicit-function**
  derivative `Leaf::dprofit_droot_collar_psi` (`leaf_model.cpp:890-922`). It
  differentiates `assim_colimited_ad(ci, vcmax_, electron_transport_, gamma_,
  km_, R_d_, ...)` treating `vcmax_/electron_transport_/gamma_/km_/R_d_` and
  `gc_const` (which contains `atm_vpd_`, line 916) as **constants** — valid only
  because temperature and VPD do not depend on the operating point today.

## Answer 1 — Does PM change the complexity of the photosynthesis modules?

**The Farquhar photosynthesis equations themselves do not change at all.** They
already take leaf temperature as input and already respond to it through the
Arrhenius block. No new terms in `assim_rubisco_limited`,
`assim_electron_limited`, `assim_colimited`, or `electron_transport`.

The complexity lands *around* photosynthesis, in three places:

1. **`Tleaf` becomes an output, not an input — a new but *non-iterative*
   coupling.** Each candidate operating point has a hydraulically-pinned `E`,
   hence an explicit `Tleaf = Tair + (Rn − λE)·ra/ρcp`, hence its own
   `vcmax_/jmax_/gamma_/kc_/ko_/R_d_`. Because `A(Tleaf)` does not feed back into
   `E`, this is a single forward pass — no iteration — provided longwave `Rn` is
   linearised or fixed (its only nonlinearity is `Rn(Tleaf)` via Stefan-Boltzmann,
   and `es(Tleaf)` if a leaf-to-air VPD is used). The real cost is that the
   "compute temperature params once and cache" optimisation (`photo_temp_cached_`,
   `leaf_model.cpp:244-253`) is **defeated**: since `Tleaf` varies per candidate
   psi, the Arrhenius block must move from once-per-`set_physiology` into the
   per-operating-point evaluation, running on every golden-section iteration of
   `find_root_collar_psi` (and every TF24f profit eval). The added work is small
   in absolute terms — exactly **7 `exp` calls** (`leaf_model.cpp:246-252`:
   `vcmax_`/`jmax_` peaked-Arrhenius = 2 each, `gamma_`/`ko_`/`kc_` = 1 each;
   `R_d_`/`km_` are arithmetic; `electron_transport_` is already recomputed every
   call). Because `Tleaf = f(E)` and `E = transpiration(psi_stem, psi_upstream)`
   is the hydraulic supply (independent of `ci`), the block is computed **once per
   operating point, before the `ci` root-find — not inside it** (mis-placing it
   inside the `psi_stem_to_ci` solve, up to `ci_niter=1000` iterations, would be
   the only way to make it expensive). A single profit evaluation already runs
   several nested root-finds (`find_psi_stem_from_psi_root`, `psi_stem_to_ci`) and
   an `E_from_Soil_to_Root_Collar` per-soil-layer loop with spline evals, which
   dominate; 7 transcendentals on top is **likely low single-digit % per solve —
   to be confirmed with the `profile-plant` skill**, not the ~10-30% I first
   guessed. `set_physiology` itself is dominated by one-time spline *setup*
   (`setup_transpiration`/`setup_root_vulnerability` over ~100 knots), which does
   **not** move and stays once-per-light-point. No caching is recoverable (`Tleaf`
   continuous in psi).

2. **The E↔gs closure.** `stom_cond_CO2` (`leaf_model.cpp:1166`) inverts
   `E = gs·VPD`. In the minimal cut this is left as-is (keep prescribed
   `atm_vpd`); only the full version replaces the VPD term with a leaf-to-air
   gradient `es(Tleaf) − ea`. Either way `E` itself stays the hydraulic supply —
   PM is used only for `Tleaf`, not to re-predict `E`.

3. **TF24f's exact analytic gradient becomes materially more complex — this is
   the single biggest cost.** `dprofit_droot_collar_psi` (`leaf_model.cpp:890`)
   presently treats the kinetic params and `gc_const` as constants (passed as
   plain doubles into `assim_colimited_ad`, line 905). Under PM they all become
   functions of the operating point via `Tleaf(E(psi))`. The exact gradient must
   gain a `dTleaf/dpsi = −(λ/ρcp)·ra·dE/dpsi` branch (`dE/dpsi` is already present
   in this function as the `transpiration_from_psi.deriv` terms) propagating
   through the *entire* peaked-Arrhenius block into `A'` (and, in the full VPD
   version, through `gc`). That means making `arrh_curve`, `peak_arrh_curve`,
   `electron_transport`, and `assim_colimited_ad` AD-traceable in `Tleaf` — the
   highest-risk, highest-effort piece. **The pragmatic answer is to defer it:**
   run TF24f on the already-present finite-difference gradient
   (`use_ad_gradient=false`, `tf24f_strategy.cpp:67-101`). The FD path
   re-evaluates `profit_at_collar_psi`, which already routes through the new
   `Tleaf` physics, so it needs **zero gradient changes** and is correct by
   construction (at the ~29% cost the code comments already record for FD).

**Verdict:** photosynthesis biochemistry — unchanged. The `Tleaf` feedback is a
single non-iterative forward pass (because `E` is hydraulically pinned), so the
leaf-energy-balance addition is low-moderate. Defeating the Arrhenius cache is a
trivial code change adding 7 transcendentals per profit evaluation — a likely
low-single-digit-% per-solve cost (confirm with `profile-plant`), dominated by the
existing hydraulic root-finds, with no correctness risk. Full drivers/traits +
RcppR6 regen — moderate, mechanical.
**TF24f's exact analytic gradient is the one high-risk/high-effort piece and
should be deferred behind the finite-difference fallback.** Complexity is
concentrated in the leaf *energy balance* and *optimiser coupling*, not in the
photosynthesis functions.

## Answer 2 — Additional environmental variables needed

Current TF24 drivers: `PPFD`, `rainfall`, `atm_vpd`, `ca`, `leaf_temp`,
`atm_o2_kpa`, `atm_kpa`, plus the multi-layer soil water state → `psi_soil`.

**Genuinely new environmental driver: only wind speed.**

| Need | Source | New? |
|---|---|---|
| Air temperature `Tair` | Reinterpret the existing `leaf_temp` driver as `Tair`; `Tleaf` becomes computed | No new driver (semantic change) |
| Above-canopy wind speed `U₀` | New extrinsic driver on `TF24_Environment` (or fixed param initially) | **Yes — the only new driver** |
| Net radiation `Rn` (shortwave) | Derived from existing absorbed PAR (`radiation`, `tf24_strategy.cpp:447`); `SW_abs ≈ 2·PAR_abs` | No |
| Longwave `Rn` component | From `Tair` (have it) + atmospheric emissivity; or fixed −40 W m⁻² offset | No |
| Actual vapour pressure / RH | Derived: `ea = es(Tair) − VPD` from existing `atm_vpd` + `Tair` | No |
| Atmospheric pressure `P` | Existing `atm_kpa` driver (scales `γ`) | No |

**New parameters (not drivers):**
- `TF24_Pars` trait: leaf characteristic dimension `d` (for `ra`); optionally
  leaf shortwave albedo `α_sw`.
- `TF24_Environment` / canopy: wind extinction coefficient `α_w`; clear-sky
  atmospheric emissivity `εa` (or derive from humidity).
- Hardcoded physical constants (like the existing block in `leaf_model.h`):
  latent heat `λ`, psychrometric constant `γ`, Stefan-Boltzmann `σ`, leaf
  emissivity `εl`, volumetric heat capacity `ρcp`.

Any interface change (new driver, new `TF24_Pars` field) must go through the
RcppR6 flow: edit the C++ header + `inst/RcppR6_classes.yml`, then
`make RcppR6 && make full_compile && make roxygen`. Never hand-edit generated
files (`R/RcppR6.R`, `src/RcppR6.cpp`, etc.).

## Answer 3 — The optimisation stays "choose the optimal root-collar ψ"

**PM does not change the structure or dimensionality of the leaf optimisation.**
Because `E` is hydraulically pinned at each operating point and `Tleaf` is an
*explicit* function of `E` (no PM inversion, and `A(Tleaf)` does not feed back
into `E`, so no inner fixed point), the decision variable remains the single
scalar root-collar potential ψ — `find_root_collar_psi` (golden-section, TF24) or
the tracked `opt_root_psi_state` ODE (gradient ascent, TF24f). PM only changes
what happens *inside* one profit evaluation:

```
candidate ψ → E = transpiration(ψ)          [hydraulic supply — unchanged]
            → Tleaf = Tair + (Rn − λE)·ra/ρcp   [NEW, explicit in E]
            → Farquhar params at Tleaf         [NEW]
            → A, then profit = A − hydraulic_cost
```

`Tleaf` is a deterministic waypoint along the ψ→profit evaluation, not a new
unknown. `profit(ψ)` stays a 1-D objective; the optimiser and its dimensionality
are identical. What shifts is the *shape* of the profit landscape (A now bends
with `Tleaf(ψ)`), so the optimum ψ\* moves — but it is still just the optimal
root-collar ψ. This confirms the issue's "no architectural change to the
optimisation structure" for TF24 specifically.

## Backward compatibility — the non-PM model as a subset

The current (non-PM) behaviour is contained two complementary ways:

1. **Runtime configuration (the mechanism to ship).** The `use_energy_balance_`
   gate (Step 1), **default off**: off → today's path runs (prescribed
   `leaf_temp`, single-shot Arrhenius, `photo_temp_cached_` intact); on → PM
   computes `Tleaf` per operating point. With the flag off, `scientific_version`
   stays put (no forced logpile cache invalidation) and PM is opt-in, matching the
   staged-hypothesis-test framing.

   **Backward-compatibility is not automatic — it splits by cut and must be
   verified:**
   - *Minimal cut (Steps 1-3)* — adds only the internal flag and internal methods;
     the Arrhenius extraction is a **pure refactor** (same arithmetic, same order),
     so with PM off the output is intended to be **bit-identical** and all existing
     runs/tests should pass unchanged. This is a design goal to **verify by running
     the full suite and diffing the saved baselines** (`results_high.RDS`,
     `tests/testthat/test_data/scenario_baseline.rds`), not an assumption — the
     refactor must preserve operation order exactly. Keep the toggle a plain C++
     `Leaf` member (default `false`), **not** a `Control` field, so the minimal cut
     touches nothing R-exposed.
   - *Full cut (Step 4)* — adding the `d` `TF24_Pars` trait and the `wind`
     extrinsic driver changes the **R-facing interface** (serialised strategy,
     `TF24_Parameters()` fields). Runtime behaviour with PM off is unchanged, but
     any snapshot/interface/round-trip test that encodes the parameter set will
     change and its reference values must be regenerated. "All tests pass
     unmodified" is therefore true only for the minimal cut.

2. **Mathematical limit (the physical interpretation).** PM's `Tleaf` reduces to
   `Tleaf = Tair` as **`ra → 0`** (boundary-layer conductance `g_Ha → ∞`:
   perfect leaf–air coupling). So the status-quo assumption is the strong-coupling
   limit of PM — useful for the sensitivity gate (sweep coupling down from ∞ and
   watch the `Tleaf − Tair` deviation appear continuously), but it will not
   reproduce a *prescribed* `leaf_temp ≠ Tair` and is not bit-identical. Document
   it as interpretation, not as the compatibility path.

## Reconciliation with the ATLS leaf thermal model

Source: `~/Downloads/ATLS_leaf_demo` (theory in `theory.qmd`, C++ in
`src/leaf_thermal_system.hpp`). ATLS is a fast (diurnal-timescale) ODE model of
leaf thermal **damage / repair / acclimation**, with states `T_L, D, I, A`
(`N = 1 − D − I`). PM and ATLS overlap on `Tleaf` and ATLS adds a layer PM lacks;
they are complementary, not competing:

1. **At the energy-balance level, ATLS's "detailed" model *is* PM.** The
   alternative balance in `theory.qmd` (attributed to Sperry et al. 2017),
   `T_L = T_A + (R_abs − εσT_A⁴ − λE)/(C_p(g_r + g_Ha))`, is algebraically the
   same steady-state balance as PM's `Tleaf`. Differences are convention only:
   conductance (`g_r + g_Ha`) vs resistance (`ra`); longwave written explicitly
   vs lumped into `Rn`; molar `C_p` vs volumetric `ρcp`. The aerodynamic submodel
   even matches — ATLS's `g_Ha = 0.189·(u/d)^0.5` is the inverse of PM's
   `ra = C·√(d/u)`, with the **same** characteristic dimension `d`.

2. **What ATLS runs in code is the *dynamic* form; PM is its steady state.** The
   C++ (`leaf_thermal_system.hpp:134`) uses a relaxation ODE
   `dT_L/dt = k_H(T_air − T_L) − g_max·S_tr(T_L)` carrying thermal inertia. PM's
   algebraic `Tleaf` is its `dT_L/dt = 0` equilibrium. Leaf thermal time constants
   (seconds–minutes) ≪ TF24's demographic step, so the quasi-steady PM `Tleaf` is
   the correct reduction; the ODE only matters if TF24 resolves sub-daily
   temperature (it does not today).

3. **ATLS's novel content is orthogonal to PM and sits downstream of it.** The
   protein damage/repair/acclimation states feed back to photosynthesis by scaling
   electron transport: `jmax_N = vcmax·J_V_ratio·N` (`leaf_thermal_system.hpp:102`).
   PM has none of this. Division of labour: **PM (= ATLS detailed balance) supplies
   `Tleaf`; ATLS consumes `Tleaf` and returns a damage multiplier `N` on
   `jmax`/quantum yield.** TF24 and ATLS already share near-identical Farquhar code
   (colimited `Ac`/`Aj`, non-rectangular-hyperbola electron transport), so the
   coupling point is unambiguous: TF24's `electron_transport()` / `jmax_`. PM feeds
   the *reversible* Arrhenius response of `jmax_`; ATLS adds a second,
   hysteretic/*irreversible* multiplier `N(Tleaf history)` on the same `jmax_`.

**Catch is architectural, not physical.** ATLS is a fast 4-state ODE resolving
diurnal heatwaves; TF24 is quasi-steady at demographic timescales with constant
`leaf_temp`. Integrating ATLS damage into TF24 needs either (a) sub-daily
temperature forcing + multiscale integration of the fast damage states within each
slow step, or (b) a steady-state/averaged damage response as a function of the
thermal regime. Either is a larger step than PM — and **PM is the prerequisite for
both**, since it makes a physically-consistent `Tleaf` available to drive the
damage kinetics.

### Nested model hierarchy (design frame)

| Model | `Tleaf` treatment | Extra states |
|---|---|---|
| Current TF24 | `Tleaf = Tair` (the `ra → 0` limit) | none |
| TF24 + PM | algebraic energy balance (= ATLS "detailed") | none |
| TF24 + PM + ATLS | same balance, optionally dynamic | protein damage `N/D/I` + acclimation `A` |

Each row is a strict reduction of the one below it. The `use_energy_balance_`
gate selects row 1 vs row 2; an eventual `use_thermal_damage_` gate + damage
states would select row 3.

## Recommended implementation (minimal core first)

Follow the doc's Section 8 "minimal first implementation" and Section 7 staged
evaluation. Both TF24 and TF24f inherit the change through the shared `Leaf`
submodel and the single `set_physiology` call site.

### Step 1 — Energy-balance functions + refactor the Arrhenius block
Add to `leaf_model.h`/`.cpp`:
- physical-constant block (alongside `leaf_model.h:20-80`): `λ`, `ρcp`, `σ`, `γ`,
  Tetens `es`/`Δ` coefficients;
- `double Leaf::leaf_temp_from_E(double E) const` → `Tair_ + (Rn_ − λ·E)·ra_/ρcp`
  (explicit; first cut uses fixed/linearised longwave in `Rn_`);
- **extract** the inline Arrhenius block (`leaf_model.cpp:246-258`) into
  `void Leaf::update_temperature_dependent_params(double Tleaf)` setting
  `vcmax_, jmax_, gamma_, kc_, ko_, R_d_, km_, electron_transport_`;
- a gate `bool use_energy_balance_` (Leaf member / Control flag) so the non-PM
  path (and Sperry/Medlyn) keeps today's single-shot behaviour.
Unit-test `leaf_temp_from_E` and the `es`/`Δ` helpers against hand computation.

### Step 2 — Wire `Tleaf` into the per-operating-point solve
- In `set_physiology` (`leaf_model.cpp:175`): reinterpret the incoming
  `leaf_temp` arg as `Tair_`; compute `Rn_` once from `radiation`
  (`tf24_strategy.cpp:447-450`), using `Tleaf=Tair` for the small LW term
  (doc Section 3.3); set a **fixed `ra_`** (doc fallback 50 s m⁻¹).
- Inside `set_leaf_states_rates_from_psi_stem` (`leaf_model.cpp:1266`) — the
  natural per-operating-point entry — on the PM path do:
  `E = transpiration(...)` → `Tleaf = leaf_temp_from_E(E)` →
  `update_temperature_dependent_params(Tleaf)` → `psi_stem_to_ci(...)` →
  `assim_colimited_`. This is the cache-bypass; `photo_temp_cached_` is skipped
  when `use_energy_balance_`.
- Keep the prescribed `atm_vpd` in `stom_cond_CO2` for the minimal cut (do not yet
  make VPD depend on `Tleaf`).

### Step 3 — TF24f gradient (defer the analytic extension)
Run TF24f with `use_ad_gradient = false` (`tf24f_strategy.cpp:67-101`); the FD
path re-evaluates `profit_at_collar_psi`, which already routes through the new
`Tleaf` physics, so **no gradient code changes are needed** and it is correct by
construction. Extend `dprofit_droot_collar_psi` with the `dTleaf/dpsi` chain
(Answer 1, point 3) only as a later, separate task if FD cost proves
unacceptable.

### Step 4 — New driver/params + regenerate (full version only)
Add `wind_speed` (`U₀`) as an extrinsic driver on `TF24_Environment`
(`tf24_environment.h:82-88` + `get_*` wrapper) and `d` (leaf dimension) as a
`TF24_Pars` field, both via `inst/RcppR6_classes.yml`, then regenerate
(`make RcppR6 && make full_compile && make roxygen`). Compute `ra = C·√(d/U₀)` in
`set_physiology` and wire the new Leaf ctor args from pars in `prepare_strategy`
(`tf24_strategy.cpp:814-819`). Only needed once the fixed-`ra` minimal core is
validated. The minimal cut in Steps 1-3 touches only `leaf_model.cpp/.h` and the
call site — **no generated-code regeneration**.

### Step 5 — Evaluation protocol (doc Section 7.5) — gate before keeping anything
Factorial Fick-vs-PM across `Tair`×`VPD`×`PAR`; report `ΔTleaf, ΔE, ΔA, Δgs`.
Keep the PM core only if `Tleaf − Tair > 2 °C` and annual carbon gain shifts
>5 % under representative conditions. Test the wind profile and the LW term
separately; drop each that fails to move outputs materially. Use the existing
`notes/plan-tf24-scenario-framework.md` scenario harness / scorecard as the
comparison infrastructure.

### Step 6 — Leaf-level demo vignette (deliverable)
A new document demonstrating leaf-level behaviour, contrasting **non-PM (default)
vs PM**. Lives in a **staging folder `overstorey_staging/`** at the package root
(e.g. `overstorey_staging/TF24_leaf_PM_demo.qmd`) — a work-in-progress area, kept
out of the built package/vignette tree for now. Can be promoted to
`inst/reports/` + a `pkgdown/_pkgdown.yml` article (cf.
`inst/reports/FF16_report.Rmd`) once the science is settled.

- **Driver:** use the already-R-exposed `Leaf` class (`inst/RcppR6_classes.yml`;
  `set_physiology(...)`, field accessors, solve entry point) — same pattern as
  `tests/testthat/test-leaf.r`. Drive a single leaf across an environment grid and
  read `profit_`, `opt_psi_stem_`, `assim_colimited_`, `stom_cond_CO2_`,
  `transpiration_`, and (PM on) `Tleaf`.
- **Interface dependency:** the demo requires `use_energy_balance_` (and the PM
  inputs `Rn`/`ra`/wind) to be **R-settable** — so it cannot run against the
  purely-internal minimal-cut flag. Expose the flag as an RcppR6 field; pair this
  step with Step 4's regeneration. (This is the one reason the demo needs the flag
  on the R interface, unlike the internal-only minimal cut.)
- **Content:**
  1. *Profit anatomy* at a single environment — assimilation benefit, hydraulic
     cost, and net profit as functions of the **root-collar potential**
     (Answer 3's decision variable, not `psi_stem` — `psi_stem` is *derived*
     from each candidate collar potential via `find_psi_stem_from_psi_root`,
     the same transport the solver itself uses), non-PM vs PM overlaid, with
     the solver's actual optimum marked on each via `evaluate_root_collar_psi`
     (so it lands on the visible curve rather than a `psi_stem` scan that
     assumes `psi_upstream == psi_soil`, i.e. negligible root resistance).
  2. *Optimal outcomes across environments* — sweep `Tair`, `VPD`, `PAR` (and
     optionally `psi_soil`); for each, plot `Tleaf − Tair`, `opt_psi_stem`, `gs`,
     `E`, `A`, and `profit`, non-PM vs PM.
  3. *Where PM matters* — highlight the hot / high-radiation / high-VPD corner,
     tying the figure back to the Step 5 evaluation gate.
- **Reuse:** `ggplot2` + `patchwork` plotting idioms from
  `R/TF24_plot_diagnostics.R`; the leaf-driving setup from `test-leaf.r`.

### Step 7 — (Deferred, separate epic traitecoevo/plant#566) ATLS thermal-damage layer
Not part of PM. Once PM supplies `Tleaf`, an optional `use_thermal_damage_` path
adds the ATLS protein damage/repair/acclimation states and multiplies `jmax_` /
quantum yield by the functional fraction `N` at the coupling point
(`electron_transport()`). Requires resolving the timescale mismatch first
(sub-daily temperature forcing + multiscale integration, or a steady-state/
averaged damage response) — a larger architectural step gated behind PM. Tracked
here only to record that PM is its prerequisite (see reconciliation section).

## Proposed new epic — ATLS leaf thermal-damage model (separate from PM)

Posted as **traitecoevo/plant#566** (`epic` label), companion to the PM epic #523.
Draft body retained below for reference:

---

**Title:** [TF24 hydraulics] Leaf thermal damage, repair & acclimation (ATLS)

**Summary.** Add a mechanistic leaf **thermal damage / repair / acclimation**
layer (the ATLS model) on top of the Penman-Monteith leaf temperature from #523.
PM supplies a physically-consistent `Tleaf`; this epic consumes `Tleaf` to track
protein-pool damage and downscale photosynthetic capacity, capturing
heat-driven loss of function and acclimation that reversible Arrhenius scaling
cannot represent.

**Depends on:** #523 (Penman-Monteith). PM is a hard prerequisite — the damage
kinetics are driven by `Tleaf`, which only becomes physically meaningful once the
energy balance is in place.

**Model (from `notes/penman-monteith/` reconciliation / `ATLS_leaf_demo`).**
Lumry-Eyring protein scheme with states functional `N`, reversibly damaged `D`,
irreversibly damaged `I` (`N+D+I=1`) plus an acclimation state `A` that raises the
critical temperature `T_crit(A)`. Temperature- and acclimation-dependent
Arrhenius rates with smooth damage/repair switches. Coupling to carbon:
photosystem-II damage downscales electron transport via `jmax → jmax·N` at TF24's
`electron_transport()` / `jmax_`. PM feeds the *reversible* Arrhenius response;
ATLS adds the *hysteretic/irreversible* `N` multiplier on the same term.

**Key challenge — timescale mismatch.** ATLS is a fast (sub-hourly) ODE resolving
diurnal heatwaves; TF24 runs quasi-steady at demographic timescales with constant
`leaf_temp`. Reconciling requires one of:
- (a) sub-daily temperature forcing + multiscale/operator-split integration of the
  fast damage states within each slow demographic step; or
- (b) a steady-state / time-averaged damage response `N(thermal regime)` applied
  as a downscaling factor (cheaper, loses hysteresis detail).

**Proposed tasks**
- [ ] Port the ATLS `Leaf` thermal system (states `T_L, D, I, A`) into the
      `plant`/odelia ODE framework as a standalone, testable submodel.
- [ ] Replace ATLS's demo relaxation-ODE `Tleaf` with the PM/Sperry energy
      balance from #523 (they are the same balance in conductance form).
- [ ] Decide and implement the timescale-reconciliation strategy (a) or (b).
- [ ] Wire the `N` multiplier into TF24 `electron_transport()`/`jmax_` behind a
      `use_thermal_damage_` gate (default off; backward-compatible like #523).
- [ ] Add damage/acclimation traits to `TF24_Pars` (RcppR6 regen) and any new
      driver (diurnal temperature amplitude) to `TF24_Environment`.
- [ ] Evaluation gate: does the damage feedback materially change annual carbon
      gain / competitive outcomes under representative Australian heatwaves,
      vs PM-only? Keep only if material (mirror #523 Section 7.5 protocol).

**Nested-model framing.** Current TF24 (`Tleaf=Tair`) ⊂ TF24+PM (#523) ⊂
TF24+PM+ATLS (this epic) — each a strict extension of the previous.

---

## Critical files
- `inst/include/plant/leaf_model.h` / `src/leaf_model.cpp` — PM functions,
  energy balance, `set_physiology`, `stom_cond_CO2` (1164), Arrhenius cache
  (244-253), `dprofit_droot_collar_psi` (890).
- `src/tf24_strategy.cpp:440` — the single `set_physiology` call site (pass
  `Tair`, wind, `Rn` inputs).
- `src/tf24f_strategy.cpp:50` — `solve_leaf`; gradient path selection.
- `inst/include/plant/models/tf24_environment.h:82-88` — driver registration
  (add `wind`).
- `inst/include/plant/models/tf24_strategy.h` (`TF24_Pars`) +
  `inst/RcppR6_classes.yml` — new `d` trait.

## Tests (sufficient coverage)

Extend `tests/testthat/test-leaf.r` (leaf-level) and the TF24/TF24f strategy and
scenario tests. Coverage must span four categories:

1. **Unit — new energy-balance maths (exact).** `es(T)`, `Δ(T)` (Tetens),
   `net_radiation_leaf(...)`, and `leaf_temp_from_E(E)` against hand-computed
   values; and the invariants: `Tleaf → Tair` as `ra → 0`, `Tleaf < Tair` under
   positive transpiration, `Tleaf` monotone decreasing in `E`.
2. **Regression / backward-compat (the key guarantee).** With
   `use_energy_balance_` **off**, assert leaf outputs (`profit_`, `opt_psi_stem_`,
   `assim_colimited_`, `stom_cond_CO2_`, `transpiration_`) are **bit-identical**
   to a recorded pre-change snapshot, and that a full SCM run reproduces
   `results_high.RDS` and `tests/testthat/test_data/scenario_baseline.rds`. This
   is what backs the "all existing tests pass" claim — it must be an actual test,
   not an assumption.
3. **PM-on behaviour.** Across an environment grid (`Tair × VPD × PAR`, plus a
   couple of `psi_soil` values): the leaf solve still returns a finite optimum
   (`opt_psi_stem_` finite, `profit_` finite); `Tleaf − Tair` grows with PAR and
   with lower `ra`; `A`/`gs` respond in the expected direction; and PM vs non-PM
   differ only where expected (≈0 in cool/low-radiation cells, material in the
   hot/high-VPD/high-PAR corner). Include a TF24f test that the FD-gradient path
   runs with PM on and its tracked `opt_root_psi_state` converges to the TF24
   golden-section optimum within tolerance.
4. **Numerical guards.** PM path stays finite under the drought / residual-soil
   edge cases that motivated #549 (non-finite `psi_soil`) — reuse the existing
   drought scenario rows so PM does not reintroduce blow-ups.
5. **Demo smoke test.** A fast test that the `overstorey_staging/` demo's
   leaf-driving code executes end-to-end on a small grid (guarded/`skip_on_cran`;
   assert it produces finite outputs for both PM on and off), so the vignette can
   never silently rot.

## Verification
1. `make RcppR6 && make full_compile && make roxygen` clean; new `TF24_Pars`
   field, `wind` driver, and the R-settable `use_energy_balance_` flag visible
   from R.
2. `devtools::test()` green, including the categories above; the backward-compat
   regression tests (category 2) pass with PM off.
3. `PLANT_RUN_SCENARIOS=1 devtools::test(filter="scenario")` still green; record a
   PM-on scorecard and diff against the pre-PM baseline
   (`notes/plan-tf24-scenario-framework.md`).
4. Sensitivity factorial (Step 5) run and effect sizes reported, and the
   `overstorey_staging/` demo renders, before any component is declared an
   improvement.
