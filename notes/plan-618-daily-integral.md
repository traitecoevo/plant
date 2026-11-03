# TF24: one-point (CODA) vs diurnal daily integral — speed and outcome test

**Experiment branch `experiment/618-daily-integral` — not for merge.** Tracks #618. The analysis is `notes/618-daily-integral/analysis.qmd`; copy it out before anything from here merges.

## Status

- [x] §1–§4 plant: Control option, environment geometry, daily integral, refusals
- [x] Verification 1–2: default bit-identical to develop (full state history, soil, offspring production; bare and seasonal stands); new tests pass; model-version snapshot gains only the four Control keys (no `scientific_version` bump, as #661 did for `absorb_epsilon`)
- [x] §5–§7 analysis stage 1: sites, c̄, 10-year runs, speed (`notes/618-daily-integral/analysis.qmd`)
- [x] Light-dependent one-point c (`one_point_c_*_slope`, `one_point_c_light_min/max`) and a 5-point Gauss–Lobatto method (`daily_integration = "lobatto"`), added 2026-11-03 at the user's request; the Lobatto idea is from `Falster-stomatal_analytical_analysis/notes/fast daily integration at fixed lambda.md`
- [x] §8 stage 2: 50-year meta-population patches and resident fitness (invader landscapes not yet run)

### Decisions taken during implementation

- **Declination is `asin(0.39785 sin(...))`.** The expression is sin(δ) (0.39785 = sin 23.44°). plant's R `solar_angle()` uses it as δ itself, ~0.6° low at the solstices; left alone there, fixed in `tf24_diurnal.h`.
- **The `DailyIntegration` enum is at namespace scope**, so a rebind between scalars copies one type.
- **The crown light is computed once per evaluation**; only PPFD_above and VPD move through the day, so the diurnal arm pays for leaf solves, not repeated crown integrals.
- **The #516 marginal-leaf gate needs no refusal**: it reads the daily `assimilation_`.
- **Stand:** lma 0.1978791 (TF24 hyperpar centre), seed rain 20 m⁻² yr⁻¹, `patch_type = "fixed"`. **Cohort schedule = plant's default ∪ a monthly grid**, the same for every arm: `run_to()` stops only at introductions, and the default schedule is years apart late in a 10-year patch.
- **Both synthetic sites persist under F**, so no rainfall adjustment was needed (temperate 0.9, arid 0.3 m yr⁻¹).

### Findings with all six methods (2026-11-03, supersede the stage-1 list below where they differ)

- **Gauss–Lobatto 5 (3 solves, no calibration) is the method to use**: leaf RMS 1.4–2.8%, 10-year stands within ~2% of D, 50-year temperate fitness +6%; 3.4–3.8× the one-point cost, 1/15 of D's.
- **Light-dependent c**: c(I) = a + b log(I/500) of intercepted noon PPFD, the same at both sites (b = 0.158). Leaf RMS 4.5–5.5%, stand totals within 2–7%, but 31–71% more plants above 1 m and 50-year fitness −66%: the log-linear form over-credits I < ~80 µmol/m2/s.
- **Constant c depends on its training sample**: widening the sample to 5% openness moved c̄ from 0.72 to 0.68 (temperate) and turned −3% biomass into −25%. 50-year fitness ~0 against D's 0.015.
- Today's TF24 (F): 50-year fitness 180× D's.
- The arid stand never reaches reproductive size in 50 years (fitness ~1e-14 in every arm), so fitness is compared at the temperate site only.

### Stage 1 findings (2026-11-03)

- c̄ (carbon, water) = 0.72, 0.74 temperate; 0.78, 0.79 arid. Noon × photoperiod overstates the day by 25–40%.
- **c depends on light**: 0.58 at canopy openness 0.2, 0.80 at full light (arid 0.64 → 0.88); season and height matter much less. A constant c over-credits shaded leaves.
- At year 10 Oc matches D on above-ground mass (−3% temperate, +7% arid), but LAI is +21% and +37%: the understorey keeps leaf area that the diurnal integral cannot pay for. O1 and F nearly double temperate biomass.
- The diurnal arm costs 46× Oc (temperate, 141 s vs 3.1 s) and 53× (arid, 118 s vs 2.2 s), tracking the leaf-solve count (49×, 56×): about 4.3 µs per solve in every arm. A first arid D run took 1141 s under machine contention; the clean rerun is the number.

## Context

TF24 turns one leaf solve into an annual carbon and water rate by multiplying by a hard-coded `60*60*12*365` (a 12 h day, every day): `tf24_strategy.h` L3160 (assimilation), L2435 (soil consumption), L2620 (sapwood controller), with a bare `//TODO: one point constant ratio and integral width for daylength` at L3072. The drivers are read once, at ODE time t; `SWdown`'s default (900 W m⁻²) is documented as a clear-sky *midday* figure. So TF24 already runs a one-point method, but with no daylength and no calibration.

The Towers et al. manuscript (`Towers-onepoint/.../main.tex`) argues that a solar-noon solve × photoperiod × a constant c̄ recovers the 15-min daily integral of profit (r² > 0.95; calibrated bias ≈ 5%), and claims this is "unlikely to significantly impact emergent ecological dynamics". plant issue #618 (diurnal closure, open) asks to "quantify the Jensen error before building anything". This test does both: same stand, same forcing, two integration methods, compare cost and emergent dynamics.

Note for the paper: Box 1 defines c = P⁽ⁿ⁾/P\* (Eq. 3, and Fig. 3's caption agrees) but then corrects with P ≈ c·P⁽ⁿ⁾ (Eq. 4), which is the inverse; it also says "mean c̄" where Fig. 3 says "median c". Here c̄ is defined directly as the multiplier, c̄ = ΣP\*/ΣP⁽¹⁾.

## Decisions (user, 2026-11-03)

- **Method is chosen by a `Control` option.**
- **One-point arm = midday amount × a constant c** (no linear/covariate calibration). Separate constants for carbon (c_A) and water (c_E), as the paper calibrates P and T separately.
- **Sites: two synthetic test cases now, temperate and arid.** Real flux-site drivers will come from Isaac later, so the driver interface must accept per-day values fitted from half-hourly data.
- **plant work stays on an experiment branch** (no PR, no merge).
- **The plan and the analysis both live on that branch.** The analysis is one self-contained qmd, to be copied out elsewhere before anything merges.

## Design — plant (experiment branch)

### 1. Control option

Add to `inst/include/plant/control.h`, defaults in `src/control.cpp`, entries in `inst/RcppR6_classes.yml` (`Control: list:`), regenerate with `make RcppR6`, add keys to `tests/testthat/test-control.R` "Defaults". Pattern: the `absorb_epsilon` change (39319c5e).

- `daily_integration` (string): `""` ≡ `"fixed-12h"` (today's code path, untouched), `"one-point"`, `"diurnal"`.
- `daily_integration_points` (int, default 48): midpoint-rule nodes across [sunrise, sunset] in `"diurnal"`.
- `one_point_c_assimilation`, `one_point_c_transpiration` (double, default 1.0).

Resolved once in `TF24_Strategy::prepare_strategy` next to `shading_model` (`strat.h` ~3658) and copied in the rebind constructor (~1556-1600).

### 2. Environment: site geometry and diurnal shapes

In `TF24_Environment` (`tf24_environment.h`), plus a small header `inst/include/plant/models/tf24_diurnal.h` for the pure functions:

- Field `latitude` (deg, default NaN; the non-default modes refuse while it is unset). Field `vpd_peak_hour` (solar hours, default 14.5).
- New driver `atm_vpd_diurnal_amp` (relative amplitude a ∈ [0, 1), default constant 0). It is a driver, not a field, so real data can vary it day to day.
- **Semantics: `SWdown(t)` and `atm_vpd(t)` are the values at solar noon of day t.** All three modes therefore see identical noon conditions; they differ only in how the day is integrated.
- Geometry: declination δ(t) uses the same formula as plant's R `solar_angle()` (`R/solar_model.R`), with day = 365·frac(t) and t = 0 ≡ 1 Jan, kept continuous (no `floor`) so rates stay smooth in t. Daylength = (24/π)·acos(−tanφ tanδ); refuse polar day/night.
- `SWdown_at(t,h) = SWdown(t) · max(0, sinβ(t,h)) / sinβ(t,12)`, which is the clear-sky geometric shape (the paper used bigleaf's potential radiation).
- `vpd_at(t,h) = atm_vpd(t) · [1 + a cos(2π(h − h_pk)/24)] / [1 + a cos(2π(12 − h_pk)/24)]`.

### 3. Strategy: the daily integral

In `TF24_Strategy`'s leaf-physiology block (`strat.h` ~2895-3070):

- Lift the shading-model block (crown-centre / MeanLight / DeepCrown) into a lambda `solve_crown(SWdown, vpd)`. `optimise_at` and `PPFD_above`/`NIR` currently read `environment.get_SWdown()`/`get_atm_vpd()`; they read local values instead. In `fixed-12h` those locals are the getter values, so this path is bit-identical.
- `"one-point"`: one `solve_crown` at noon. Then scale `profit_` by c_A, and `transpiration_` and every `soil_consumption_[a]` by c_E.
- `"diurnal"`: for k = 0..n−1, h_k = sunrise + (k+½)·Δ, Δ = daylength/n; call `solve_crown(SWdown_at, vpd_at)` and accumulate `profit_`, `transpiration_` (which `shadow_cost()` reads), every `soil_consumption_[a]`, and the diagnostics as daylight means (sum/n). The DeepCrown branch (3024-3062) is the template for which members to integrate. **n = 1 puts the single node at solar noon, so it must reproduce one-point with c = 1.**
- Annual conversion: in the two new modes, L3160 and L2435 use a `daylight_seconds_per_year_ = daylength_h(t)·3600·365`. **`fixed-12h` keeps the literal `60*60*12*365` chains** (the DO-NOT-HOIST note at L1905: hoisting moves offspring production by about 0.2%).
- Aux outputs in the new modes are daylight means (state this in the header comment). Dark-ish nodes near sunrise/sunset return profit ≈ −R_d, not an error (phylloptim `leaf_photosynthesis.hpp:266-311`); keep the negatives in the integral.
- phylloptim solves carry no warm start (`leaf_collar_solve.hpp:313-327`), so a sequence of solves equals fresh solves. The replay record (`leaf_solved_points`) is one list per evaluation with one entry per solve, so n solves is consistent.

### 4. Refusals (non-default modes only, in `prepare_strategy`)

Active S (AD/gradients), TF24f (single tracked collar), TF24t (`Tleaf_` drives damage, and a daylight-mean Tleaf is meaningless there), `a_sw > 0` (2 extra solves at one light), `a_pl0 > 0` if the #516 gate reads a single-light quantity (check at ~3191), latitude unset, n < 1, c ≤ 0, unknown method string. The tests cover only the essentials, since this is an experiment branch.

## Design — analysis (one qmd on the branch)

All of §5–§8 lives in `notes/618-daily-integral/analysis.qmd`: helper functions are defined in its own chunks (no package code), expensive runs are cached to `.rds` beside it (git-ignored), and figures come from `plot_*()` functions in the qmd.

### 5. Test sites (synthetic, illustrative values; not cited climatology)

An R builder `site_drivers(site, years)` produces spline drivers on 48 knots per year (as `R/scenario_eval.R` does) and sets `latitude`, `vpd_peak_hour`:

- **temperate**: lat −33.6° (Cumberland Plain, the paper's base site). SWdown_noon(t) is clear-sky geometry × transmissivity ~0.7. Noon VPD has a seasonal sinusoid with a summer peak, roughly 0.8–1.8 kPa, a ≈ 0.4. Rainfall is ~0.9 m yr⁻¹ with mild seasonality (amp_frac 0.3, summer max).
- **arid**: lat −22.3° (Alice Springs Mulga). Transmissivity ~0.8. Noon VPD is roughly 1.5–2.6 kPa, a ≈ 0.25. Rainfall is ~0.3 m yr⁻¹ with amp_frac 0.8, summer-dominant.
- **Constraint:** the leaf refuses VPD ≥ e_s(T_air) (`leaf.hpp:2401`), and `leaf_temp` is pinned at 25 °C (e_s ≈ 3.17 kPa). The builder asserts the daily VPD peak stays < 0.95·e_s(25 °C), so "arid" is as arid as the pin allows.
- **Gate:** under `fixed-12h` the default TF24 stand must persist 10 yr at both sites. If arid collapses, raise rainfall and record the value used.

### 6. Calibrating c (the paper's site-level c̄)

`calibrate_c(site)` uses TF24 individuals through plant's own code path (`env$time <- t; env$set_fixed_environment(openness); env$set_soil_water_state(θ); ind$compute_rates(env)`, as in `test-strategy-tf24.R:548-559`):

- 100 random days × heights {1, 4, 12 m} × openness {1, 0.5, 0.2} × θ at the wet/median/dry quantiles of that site's `fixed-12h` soil trajectory.
- Read assimilation and summed consumption under `diurnal` (n = 48) and `one-point` (c = 1).
- c̄_A = ΣP\*/ΣP⁽¹⁾ and c̄_E = ΣT\*/ΣT⁽¹⁾ (ratio of sums: the mean of ratios blows up where noon profit is ≈ 0 or negative, which the paper found at dry sites). Report the per-sample spread of c too.
- The same sample gives the convergence check that picks the reference n: compare n = 24/48/96 and accept 48 if the aggregate |P\*₄₈ − P\*₉₆|/P\*₉₆ < 1e-3.

### 7. Stage 1 runs: 10-year stands (per site)

`patch_type = "fixed"`, `max_patch_lifetime = 10`, one species at TF24 hyperpar defaults, MeanLight (default), spline drivers, `collect = TRUE`. Arms:

| Arm | `daily_integration` | c |
|---|---|---|
| D (reference) | diurnal, n = 48 | — |
| O1 | one-point | 1 |
| Oc | one-point | c̄_A, c̄_E (site) |
| F (status quo) | fixed-12h | — |

- **Outcomes:** LAI(t) (the `patch_lai` helper in `overstorey/guides/patch-dynamics.qmd`), stand leaf mass / above-ground mass / basal area (`stand_census`), soil water and ψ(t) (`tidy_env`), height distributions at t = 1, 2, 5, 10, growth (dH/dt vs H at t = 2, 5, 10, from `net_mass_production_dt` aux or by differencing height along nodes), cumulative stand assimilation and transpiration, offspring production.
- **Metrics vs D:** relative difference at t = 10 (LAI, biomass), max relative |ΔLAI| over t, height median and p90 at t = 10, offspring-production ratio.
- **Speed:** wall time (min of 3 interleaved reps), leaf solves (`plant:::census_operating_point_counts_tf24`), ODE steps. Report D/Oc time and solve ratios.

### 8. Stage 2: fitness (only once stage 1 shows D is affordable)

Meta-population patch at the default lifetime. Compare resident `net_reproduction_ratios` / offspring production for D vs Oc, and an invader fitness landscape over one trait (lma or stem_P50, ~7 values, `scm$run_mutant`). Does Oc preserve the landscape's shape and argmax? First confirm that `run_mutant` runs in diurnal mode.

### 9. Later: Isaac's flux-site data

Fit per-day noon SWdown, noon VPD and the relative amplitude from half-hourly data (smoothed, as the paper does), and feed them in as the same three drivers. Known gaps to resolve then:

- VPD above 3.17 kPa will occur at warm sites. That needs a `leaf_temp` driver, which reopens the 25 °C pin.
- Cloudy-day noise in SWdown.
- Day-to-day spline wiggle means more ODE steps; NEWS reports ~9–10k accepted steps on a 41-yr daily record.

## Where things live

Everything is on plant branch `experiment/618-daily-integral`, cut from `origin/develop` (bdc61011, `scientific_version = 14`):

- **Worktree:** under the session scratchpad, because sibling worktrees get pruned.
- **Plan:** `notes/plan-618-daily-integral.md`. This plan, committed first, is the single source of truth: keep it current and mark what's done.
- **Analysis:** `notes/618-daily-integral/analysis.qmd`. It records the plant commit it rendered against.
- **C++/R changes** as in §1–§4.
- **Commits** are local. **Ask before pushing** the branch to origin: pushing protects it from worktree pruning, but it is outward-facing.
- **No GitHub comment on #618** without asking.

## Build

- Scratch R library `<scratchpad>/Rlib-618` holding phylloptim **v0.26.2**: plant pins `== 0.26.2`, but 0.27.2 is installed. odelia 0.10.5 from the site library is already the pinned version.
- `make RcppR6`; `rm -f src/*.o` before each install (R CMD INSTALL silently reuses objects). Use `GIT_CONFIG_GLOBAL=/dev/null R_MAKEVARS_USER=` in the sandbox.
- The baseline arm for bit-identity is `origin/develop` built into a second scratch library.

## Verification

1. **Default untouched:** default TF24 stand (10 yr, both sites' drivers and the bare default env) on the branch vs `origin/develop`: `identical()` on the full state history and offspring production (`%.17g`). Then run `NOT_CRAN=true TESTTHAT_PARALLEL=false` `test_check("plant")`, filtered to tf24|control|strategy|scm; no snapshot moves.
2. **New `tests/testthat/test-tf24-daily-integration.R`:**
   - daylength = 12 h at δ = 0; `SWdown_at(noon) == SWdown(t)` and `vpd_at(noon) == atm_vpd(t)`; SWdown_at(sunrise) = 0;
   - **diurnal n = 1 equals one-point c = 1** (assimilation and consumption, to rounding);
   - one-point assimilation scales exactly with c_A, and consumption with c_E;
   - diurnal n = 48 vs 96 is close at a test state;
   - the solve census shows n solves vs 1;
   - refusals: latitude unset, unknown method.
3. **Analysis sanity:** F and O1 differ only by the daylength factor (check by ratio); annual mean daylength ≈ 12 h.
4. **Report** in plain language, results first: does Oc track D on LAI, size structure, growth and (stage 2) fitness, and at what speed ratio? Then ask about pushing the branch and posting to #618.
