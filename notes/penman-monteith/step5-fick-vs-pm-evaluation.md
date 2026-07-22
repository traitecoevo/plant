# Step 5 evaluation — Fick's law vs Penman-Monteith (leaf level)

Implements the doc Section 7.5.1 factorial gate: does the Penman-Monteith leaf
energy balance change leaf outputs enough, and in the predicted places, to be
worth keeping? Companion to `implementation-plan.md` (Steps 1–4 built the
feature; this is the "gate before keeping anything" measurement).

## Method

A single TF24-default `Leaf` (via `TF24_Strategy()$pars`; root defaults
`root_b/root_c/beta_R_*` from `tf24_strategy.h`) is driven across the factorial

    Tair ∈ {15, 25, 35, 40} °C  ×  VPD ∈ {0.5, 1, 2, 3} kPa
                                 ×  PAR ∈ {200, 500, 1000, 2000} µmol m⁻² s⁻¹

with Fick's law (`use_energy_balance_ = FALSE`, `Tleaf = Tair`) and with PM
(`TRUE`). For each cell we solve the root-collar profit optimum
(`find_root_collar_psi`) and record `E`, `A` (`assim_colimited_`), `gs`
(`stom_cond_CO2_`) and — on the PM path — the operating-point leaf temperature
`Tleaf = Tair + (Rn − λE)·ra/ρcp` (recovered in R from the exposed `Tair_`,
`Rn_`, `ra_`, `transpiration_`). PM inputs use the minimal-cut submodels: net
radiation `Rn = 2·PAR/4.57 − 40` W m⁻² (shortwave ≈ 2×PAR, fixed −40 clear-sky
longwave), and `ra = C_ra·√(d/U₀)` with `d = 0.05 m`, `U₀ = 2 m s⁻¹`
(⇒ ra ≈ 31.6 s m⁻¹). Prescribed `atm_vpd` is kept (leaf-to-air VPD deferred).

Leaf geometry is tuned to a well-watered, well-rooted, moderate-conductance
operating point so the profit optimum **opens** the stomata; at the shut-down
optimum the Fick-vs-PM contrast is degenerate (A ≈ 0 either way). This is a
single isolated leaf, not an SCM run — see caveats.

Regenerate: `Rscript notes/penman-monteith/step5_fick_vs_pm.R <outstem>`
(writes `<outstem>.rds` + `.csv`). Committed data: `step5_fick_vs_pm.csv`.

## Results (effect sizes)

Slice at VPD = 2 kPa (ΔTleaf in °C; ΔA, ΔE relative to Fick, %):

| PAR  | Tair | ΔTleaf | ΔA %  | ΔE %  |
|-----:|-----:|-------:|------:|------:|
|  200 |  15  | −0.06  |  +0.9 |  +2.5 |
|  500 |  15  |  +1.72 |  +1.5 |  +2.1 |
| 1000 |  15  |  +5.32 |  +7.4 |   0.0 |
| 2000 |  15  | +16.86 | −15.0 |   0.0 |
|  200 |  25  |  +0.05 |  +2.1 |  +5.7 |
| 1000 |  25  |  +5.32 | −13.4 |   0.0 |
| 2000 |  25  | +21.86 | −97.6 | −97.0 |
|  200 |  35  |  +0.42 |   0.0 |  +3.3 |
| 1000 |  35  | +10.25 | −91.5 | −90.6 |
| 2000 |  35  | +22.01 | −100  | −100  |
| 1000 |  40  | +10.42 | −95.0 | −92.0 |
| 2000 |  40  | +22.01 | −101  | −100  |

Global (all 64 cells): **max |Tleaf − Tair| = 22.0 °C**, **median 4.9 °C**;
**max |ΔA| ≈ 102 %**, **max |ΔE| ≈ 300 %**. Cool/low-light corner
(Tair ≤ 25, PAR ≤ 500): mean |ΔTleaf| = 1.15 °C, mean |ΔA| = 1.8 %.

## Verdict — KEEP the PM core

Both gate thresholds are cleared, and the effect is localised exactly where the
doc predicted:

- **|Tleaf − Tair| > 2 °C under target conditions** — yes (up to 22 °C; ≥ 5 °C
  routine once PAR ≥ 1000). ✅
- **Materially changes A** — yes; in the hot / high-radiation / high-VPD corner
  A collapses by 90–100 % as the overheated leaf shuts down, far exceeding the
  5 % bar. ✅
- **Negligible in cool/low-radiation cells** — yes (< 2 °C, < 2 % A), so PM does
  not perturb conditions where Fick's law is adequate. ✅

**Structural, not incidental.** ΔE tracks ΔgS (stomata close as the leaf
overheats and A falls); E is otherwise close to hydraulically pinned, so PM's
correction flows mainly into the **carbon** side (A via Farquhar temperature
scaling) — consistent with the doc's argument that the Tleaf→Farquhar feedback
is a structural bias correction, largest in hot/high-radiation Australian
conditions. Some moderate-temperature cells show PM *raising* A (warming toward
the Vcmax/Jmax optimum), so the correction is directional, not a uniform penalty.

## Robustness note (required by PM)

Extreme energy-balance heating raises the CO₂ compensation point (Γ*) until
`psi_stem_to_ci`'s `[Γ*, ca]` bracket contains no supply==demand root, which
previously threw. Added a **PM-only** graceful fallback in `psi_stem_to_ci`
(operate at the compensation point, gross A = 0, net A = −R_d) so the leaf shuts
down smoothly instead of crashing; gated on `use_energy_balance_`, so the non-PM
path keeps its fail-fast contract and stays bit-identical (verified).

## Caveats / what this does NOT yet establish (doc 7.5.2–7.5.4)

- **Isolated leaf, tuned geometry.** Absolute magnitudes depend on the chosen
  conductance/root mass; the *pattern* (hot-corner collapse, cool-corner
  no-op) is the robust result. An SCM-level run confirms PM is material there
  too (offspring production shifts by orders of magnitude, PM on vs off).
- **The 22 °C at PAR = 2000 is an upper bound** — radiative heating with little
  evaporative cooling, using `Rn = 2·PAR` (no NIR albedo split) and fixed
  `ra`/LW. Not yet validated against measured Tleaf.
- **Still to run before publication:** fixed-ra vs wind-profile (7.5.2), Rn ±30 %
  sensitivity (7.5.3), and the annual SCM integration (7.5.4). The wind profile
  (α_w) and the longwave term remain candidates to drop if they fail their own
  sensitivity gates.
