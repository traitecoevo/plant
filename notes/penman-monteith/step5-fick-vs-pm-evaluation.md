# Step 5 evaluation — Fick's law vs Penman-Monteith (leaf level)

Implements the doc Section 7.5.1 factorial gate: does the Penman-Monteith leaf
energy balance change leaf outputs enough, and in the predicted places, to be
worth keeping? Companion to `implementation-plan.md` (Steps 1–4 built the
feature; this is the "gate before keeping anything" measurement).

The 2026-10 re-run comes first. The July run follows unchanged, as the record of what it found. The three sensitivity tests July left open (doc 7.5.2–7.5.4) are in [`sensitivity-results.md`](sensitivity-results.md).

## 2026-10 re-run

Against plant `develop` at `bd09645d` (2.0.0.9002), phylloptim `edaffe8` (0.10.2) and odelia `88cbff7` (0.7.3), the pins in plant's DESCRIPTION.

**Verdict: the gate still passes, by a wider margin than in July. Keep the PM core, default off.** In the corner the gate names (Tair ≥ 35 °C, PAR ≥ 1000, VPD ≥ 2 kPa), the PM leaf sits 10.5–22 °C above air and shuts. Its net photosynthesis goes from +0.8 to +9.8 µmol m⁻² s⁻¹ under Fick to −3.0 to −3.8 under PM: it respires. The gate asks for 2 °C and a 5 % change in A. What changed since July is how the leaf behaves on the way there. It now transpires hard to cool itself, and when it cannot cool enough it loses more photosynthesis than it did in July. The sensitivity tests add a caution: the PM answer depends strongly on two inputs nobody has calibrated, and at TF24's default drivers it collapses the stand. So PM is worth keeping, but not ready to switch on.

### What the quantities mean

- **Tleaf − Tair**: how much warmer the PM leaf is than the air, °C. Fick holds the leaf at air temperature by construction.
- **ΔA**: the PM leaf's net photosynthesis relative to the Fick leaf in the same cell, %. Reported only where the Fick leaf's A exceeds 1 µmol m⁻² s⁻¹. Net A is gross assimilation minus dark respiration, so a leaf that shuts and respires reads below −100 %. That means "a net carbon source", not "more than total loss".
- **ΔE**: the same comparison for transpiration, kg H₂O m⁻² s⁻¹.
- **Transpiration ceiling**: the most water this leaf's roots can deliver, 6.89 × 10⁻⁵ kg m⁻² s⁻¹. That flow carries away about 169 W m⁻² as latent heat.

### Headline numbers, July against October

All 64 cells of the same factorial (Tair × VPD × PAR), same leaf, same geometry.

| quantity | July | October |
|---|---:|---:|
| max \|Tleaf − Tair\| | 22.0 °C | 22.0 °C |
| median \|Tleaf − Tair\| | 4.9 °C | 4.6 °C |
| cells where \|ΔA\| > 5 % | 41 of 64 | 43 of 59 |
| cells where PM transpires > 10 % more than Fick | 11 | 33 |
| PM cells at the transpiration ceiling | 15 | 23 |
| PM cells that shut (E < 10⁻⁶) | 12 | 16 (21 at the true optimum) |
| PM cells with net A < 0 | 4 | 20 |
| cool corner (Tair ≤ 25, PAR ≤ 500): mean / max \|ΔA\| | 1.8 % / 5.6 % | 3.8 % / 12.4 % |

The maximum warming did not move because it is the no-transpiration limit, Rn·ra/ρcp = 835 W m⁻² × 31.6 s m⁻¹ / 1200 = 22.0 °C at PAR 2000. That arithmetic is unchanged since July. Five October cells have no ΔA because the Fick leaf at 40 °C now has A ≤ 1 (point 4 below).

Slice at VPD = 2 kPa (ΔTleaf in °C; ΔA and ΔE in %; solver results, as TF24 would use them):

| PAR | Tair | ΔTleaf Jul | ΔTleaf Oct | ΔA Jul | ΔA Oct | ΔE Jul | ΔE Oct |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 200 | 15 | −0.06 | −0.11 | +0.9 | +2.3 | +2.5 | +6.1 |
| 500 | 15 | +1.72 | +1.02 | +1.5 | +4.7 | +2.1 | +25.6 |
| 1000 | 15 | +5.32 | +6.03 | +7.4 | −4.3 | 0.0 | 0.0 |
| 2000 | 15 | +16.86 | +17.56 | −15.0 | −49.6 | 0.0 | 0.0 |
| 200 | 25 | +0.05 | −0.05 | +2.1 | +5.7 | +5.7 | +14.8 |
| 500 | 25 | −0.44 | +0.26 | +0.7 | −0.9 | 0.0 | 0.0 |
| 1000 | 25 | +5.32 | +6.03 | −13.4 | −32.1 | 0.0 | 0.0 |
| 2000 | 25 | +21.86 | +21.98 | −97.5 | −121.2 | −97.0 | −99.3 |
| 200 | 35 | +0.42 | +0.46 | 0.0 | +0.3 | +3.3 | +16.8 |
| 500 | 35 | −0.44 | +0.26 | +39.4 | +32.8 | +190.2 | +187.1 |
| 1000 | 35 | +10.25 | +7.98 ⚠ | −91.5 | −147.2 ⚠ | −90.6 | +22.5 ⚠ |
| 2000 | 35 | +22.01 | +21.98 | −100.0 | −136.4 | −100.0 | −99.3 |
| 200 | 40 | +0.66 | −3.20 ⚠ | −5.5 | n/a | +1.6 | +1501 ⚠ |
| 500 | 40 | +4.43 | +2.93 | −65.7 | −98.7 | −58.8 | +345.8 |
| 1000 | 40 | +10.42 | +10.45 | −94.9 | −410.2 | −92.0 | −92.8 |
| 2000 | 40 | +22.01 | +21.98 | −101.4 | −341.5 | −100.0 | −92.8 |

⚠ marks cells where plant's solver did not return the profit maximum (point 5). At the true maximum, 1000/35 shuts at +10.5 °C with A = −3.65. At 200/40 the leaf nearly shuts at +1.1 °C with A = 0.29, not the 3.81 the solver reports.

### What moved, and why

1. **The PM leaf now transpires to cool itself.** In 33 of 64 cells it moves more than 10 % more water than the Fick leaf; July had 11. Twenty-three cells sit at the transpiration ceiling, against 15 in July. At PAR 500 and Tair 25–35 every PM leaf runs at the ceiling and stays within 0.3 °C of air; in July those leaves ran up to 3.9 °C warm. Two phylloptim changes did this. #84 put the cooling effect of transpiration into the solver's optimality condition: July's solver took the derivative at air temperature, with no dA/dTleaf term, and landed up to 0.83 MPa from the true optimum. #93 made cooling worth more: a cooler leaf also has a smaller leaf-to-air vapour deficit, so each unit of water buys more stomatal conductance.
2. **A leaf that cannot cool enough loses more photosynthesis.** At PAR 2000, Tair 15, VPD 2 the PM leaf is about 17 °C warm in both runs. Its loss in A grew from 15 % to 50 %. At 1000/25/2 it grew from 13 % to 32 %. The cause is #93: Fick's law now divides by the deficit at the leaf's temperature, 3.3–5.2 kPa in these cells, not the 2 kPa of air. July overstated stomatal conductance by that ratio wherever the leaf was warm, and overstated A with it.
3. **A shut leaf now respires.** Shut PM leaves report A = −3.0 to −3.8 µmol m⁻² s⁻¹, against about 0 in July. That is why ΔA now runs from −118 % to −410 %. #86 made dark respiration rise with temperature (Tjoelker's declining Q10, R_d_25 = 1.44). #111 and #117 made the shut-down exits evaluate respiration at the leaf's own temperature, which is 45–62 °C here, not air temperature.
4. **The Fick baseline moved as well**, independently of PM. (a) At 35–40 °C the Fick leaf's A is 0.7–1.9 µmol lower, because of the same #86 respiration: R_d is 2.59 at 35 °C and 3.17 at 40 °C. (b) The transpiration ceiling fell 14 %, from 7.99 to 6.89 × 10⁻⁵. The root network and both vulnerability curves are identical to July's (checked against `b81786c6^`). What changed is that the root collar now stops at the root's critical potential, 5.87 MPa. Today's code reproduces July's ceiling at a collar of about 6.8 MPa. *Inferred*: July's solve let the collar run past the root's critical potential. This was not traced to a single PR.
5. **The solver sometimes returns the lower of two maxima.** On the PM path, profit as a function of the collar potential can peak twice: transpire hard and cool, or shut. A brute-force scan of the collar bracket (`pm_scan_collar`, 300 points plus a golden-section refinement) beat the solver by more than 0.1 µmol m⁻² s⁻¹ in 6 of 64 PM cells and 1 of 64 Fick cells. The worst cases: at 1000/35 (all four VPDs) the solver reports profit −8.26 where shutting gives −3.86. At 200/40/VPD ≥ 2 it reports a "rescued" leaf with A = 3.7–3.8 where the optimum is nearly shut. A further 20 PM cells differ by 0.02–0.05, mostly at the zero-flow end of the bracket. None of this changes the verdict: in the gate corner net A is below zero at both the solver's answer and the true optimum. It does mean a PM result in a single cell should be checked against the scan. The re-run CSV carries both: `gap_on`/`gap_off` and the `*_global` columns.

### Reproduce

`Rscript notes/penman-monteith/step5_fick_vs_pm.R notes/penman-monteith/step5_fick_vs_pm_2026-10` from the repo root (about 1 s) writes `step5_fick_vs_pm_2026-10.{csv,rds}`; only the CSV is committed. July's output stays in `step5_fick_vs_pm.csv` for comparison. The script was ported, not redesigned. The leaf now comes from phylloptim's constructor, seated on TF24_floor at λ₀ = 0 exactly as `TF24_Strategy` seats it (identical to `find_root_collar_psi()` on the 24 cells checked). `set_physiology` takes a root network in place of root mass and no longer takes the sapwood volume. The shared driver is `pm_leaf_common.R`. Its R-side leaf temperature, `Tair + (Rn − λE)·ra/ρcp` clamped to [−40, 70] °C, matches the leaf's own `Tleaf_` (as reported through TF24's `Tleaf` auxiliary) to 1.4 × 10⁻¹⁴ °C.

## July 2026 run (superseded by the 2026-10 re-run above; kept as the record)

### Method

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

### Results (effect sizes)

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

### Verdict — KEEP the PM core

*2026-10: the first two bullets still hold. The third no longer does: the cool corner now moves A by up to 12 %. "E is otherwise close to hydraulically pinned" no longer holds either, because the PM leaf now transpires to cool. See the re-run.*

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

### Robustness note (required by PM)

*2026-10: this code now lives in phylloptim, and the fallback is no longer PM-only. Since respiration rises with temperature (phylloptim #86), the prescribed-temperature path can also reach a bracket with no root, and it now shuts down there too. Any other root-find failure still stops.*

Extreme energy-balance heating raises the CO₂ compensation point (Γ*) until
`psi_stem_to_ci`'s `[Γ*, ca]` bracket contains no supply==demand root, which
previously threw. Added a **PM-only** graceful fallback in `psi_stem_to_ci`
(operate at the compensation point, gross A = 0, net A = −R_d) so the leaf shuts
down smoothly instead of crashing; gated on `use_energy_balance_`, so the non-PM
path keeps its fail-fast contract and stays bit-identical (verified).

### Caveats / what this does NOT yet establish (doc 7.5.2–7.5.4)

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
  *2026-10: all three run; see [`sensitivity-results.md`](sensitivity-results.md).*
