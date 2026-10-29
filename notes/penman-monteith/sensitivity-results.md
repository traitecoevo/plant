# Penman-Monteith sensitivity tests (doc 7.5.2–7.5.4), 2026-10

The three tests `penman_monteith_implementation.md` §7.5 required before PM outputs are used, and which July's Step 5 left open. They were run against plant `develop` at `bd09645d` (2.0.0.9002), phylloptim `edaffe8` (0.10.2) and odelia `88cbff7` (0.7.3). The re-run of the Step 5 factorial (7.5.1) is in [`step5-fick-vs-pm-evaluation.md`](step5-fick-vs-pm-evaluation.md).

## Verdict

**PM matters far more than the #523 gate asks, and its answer depends on inputs nobody has calibrated. Keep the PM core, default off. Do not switch it on for science yet.**

- **Boundary-layer resistance (7.5.2): it matters, so a fixed value is not neutral.** Over the wind model's own range (wind 0.5–5 m s⁻¹), leaf temperature moves by several °C and photosynthesis by more than 5 % in a quarter to a half of the cells. Stand-level offspring production spans five to nine orders of magnitude between 0.8 and 5 m s⁻¹. The doc's fixed fallback of 50 s m⁻¹ is the same as a 0.8 m s⁻¹ wind. The doc's exponential within-canopy profile is not usable as written: past about 1.5 units of leaf area it leaves the doc's own 20–200 s m⁻¹ range and cools shaded leaves below the wet-bulb temperature.
- **Net radiation (7.5.3): it matters, so it cannot be simplified.** ±30 % on net radiation changes photosynthesis by more than 5 % in about half the cells. The fixed −40 W m⁻² longwave term corresponds to a humid or overcast sky. A clear sky means 76–159 W m⁻² of longwave cooling and leaves 1–3 °C cooler.
- **Stand level (7.5.4): PM is a regime change at TF24's defaults, not a correction.** PM changes one year's above-ground growth by −4 % at 10 °C and −65 % at 15 °C, and the stand dies at 30 °C. Offspring production falls 70-fold at 10 °C. Almost all of this runs through the larger vapour deficit at a warmer leaf (phylloptim #93), not through leaf temperature acting on photosynthesis. Fick at PM's leaf temperature grows normally. Fick at PM's leaf temperature *and* PM's leaf-to-air deficit reproduces the collapse.

Against the #523 gate (keep PM if Tleaf − Tair > 2 °C and A or annual carbon gain shifts by more than ~5 %): cleared at leaf level in the gate corner (10.5–22 °C, A from positive to negative) and at stand level at every air temperature from 15 °C up.

## What was run

| test | script | output | grid | runtime |
|---|---|---|---|---|
| 7.5.2 ra, factorial | `sensitivity_leaf.R` | `sensitivity_ra.csv` | 64 cells (Tair 15/25/35/40 × VPD 0.5/1/2/3 × PAR 200/500/1000/2000) × 6 ra values | ~8 s, all three leaf tests together |
| 7.5.2 ra, canopy transect | `sensitivity_leaf.R` | `sensitivity_ra_canopy.csv` | depth L = 0–2 × Tair 15/25/35 × VPD 1/2 × 3 ra treatments | |
| 7.5.3 Rn | `sensitivity_leaf.R` | `sensitivity_rn.csv` | 64 cells × (Rn ×0.7/×1/×1.3, longwave at εa 0.7/0.8/0.9) | |
| 7.5.4 stand | `sensitivity_scm.R` | `sensitivity_scm.csv` | Tair 10–35 in 5 °C steps × PM off/on, plus matched and wind runs | ~12 min on an unloaded machine, 10 of them in the one dying Tair = 35 PM-on cell (that cell took 61 min on a re-run under heavy load) |

Run each from the repo root with `Rscript notes/penman-monteith/<script>`. All of them print and check the stack (`pm_check_stack`), and all share the leaf of Step 5 (`pm_leaf_common.R`).

**Leaf runs.** Every cell has PM on and is reported at its **global** profit maximum: a brute-force scan of the root-collar bracket plus a golden-section refinement, taken wherever it beats plant's solver. The reason is that plant's solver sometimes returns the lower of two maxima: it missed by more than 0.1 µmol m⁻² s⁻¹ in 6 of 64 baseline cells (Step 5 write-up, point 5). A perturbation that only moved the solver between maxima would otherwise read as a sensitivity of the model. The `gap` column records where the scan won.

ra and Rn are injected into the leaf after `set_physiology` derives them and before the solve reads them. Both are read only by the solve, in `leaf_temp_from_E`. The script checks the injection first: injecting the model's own ra changes nothing, and injecting ra = 50 gives the same leaf, to 10⁻¹², as reaching 50 through the wind model at U0 = 0.8 m s⁻¹.

**Stand runs.** TF24's SCM, one species (lma 0.0825, hmat 5 m, birth rate 20; the 5-year fixture of `test-strategy-tf24.R`), density carried in **birth date** (`Control(node_density_in_birth_date = TRUE)`, the coordinate TF24 needs, #590). PM is switched on through each strategy's `pars$use_energy_balance`, asserted after it is set. Drivers are constant at TF24's defaults (PPFD 1800, VPD 1 kPa, wind 2 m s⁻¹, rainfall 1, ca 40 Pa), except air temperature (`leaf_temp`). For each case:

- a 1-year patch, read at t = 1: leaf area and above-ground mass per m² (TF24's `stand_census`), and the water the stand drew from the soil (`sum_resource_depletion`);
- a 5-year patch: offspring production;
- for the main grid only, a 5-year two-species patch (lma 0.0825 against 0.10): each species' offspring.

## 7.5.2 — Boundary-layer resistance

**What ra is.** The resistance of the still air around a leaf to carrying heat away, in s m⁻¹. The leaf runs Tleaf = Tair + (Rn − λE)·ra/ρcp. So ra multiplies the leaf's energy imbalance: it sets how far a leaf that absorbs more than it evaporates heats above the air. plant computes ra = 200·√(d/U0) from leaf size d (0.05 m) and above-canopy wind U0 (2 m s⁻¹ by default), which gives 31.6 s m⁻¹. There is no within-canopy profile.

**Across the factorial**, each ra value is compared with the default 31.6. "Cells that switch" counts cells that change between open and shut (shut means E < 10⁻⁶). ΔA is counted over the 36 cells where the default leaf's A exceeds 1 µmol m⁻² s⁻¹. The last column follows one representative open cell (PAR 1000, Tair 25 °C, VPD 2 kPa).

| ra (s m⁻¹) | what it represents | cells with \|ΔA\| > 5 % (of 36) | cells that switch (of 64) | median \|ΔTleaf\| (°C) | PAR 1000 / 25 °C / 2 kPa: Tleaf, A |
|---:|---|---:|---:|---:|---|
| 20 | wind model at U0 = 5 m s⁻¹ | 9 | 10 | 1.4 | 28.8 °C, 13.7 |
| 31.6 | wind model at U0 = 2 m s⁻¹ (plant's default) | — | — | — | 31.0 °C, 11.5 |
| 50 | the doc's (and the leaf's) fixed fallback; U0 = 0.8 m s⁻¹ | 13 | 9 | 3.5 | 34.5 °C, 7.8 |
| 63 | wind model at U0 = 0.5 m s⁻¹ | 17 | 11 | 6.0 | 37.1 °C, 5.1 |
| 110 | the doc's profile one unit of LAI into the canopy (α_w = 2.5) | 28 | 13 | 12.7 | 61.6 °C, shut |
| 200 | top of the doc's "realistic" 20–200 range | 31 | 17 | 17.3 | 70 °C (the clamp), shut |

So outputs are not insensitive to ra over any realistic range, and the doc's condition for dropping wind dependence is not met. Within the wind model's own range (20–63 s m⁻¹), the representative cell's photosynthesis falls from 13.7 to 5.1 µmol m⁻² s⁻¹ (−63 %). Replacing the wind model with a fixed 50 is not a neutral simplification: it asserts a still day.

**Within the canopy.** The doc's profile is U(L) = U0·exp(−α_w·L) over cumulative leaf area L, which makes ra = ra0·exp(α_w·L/2). With α_w = 2.5–3.5 that gives ra = 110–182 s m⁻¹ at L = 1, 385–1047 at L = 2 and 1345–6026 at L = 3. It leaves the doc's own 20–200 range by L ≈ 1–1.5. Along a depth transect, absorbed PAR follows TF24's light model (0.5 × 1800 × e^(−0.5L)). Three treatments were compared: today's single top-of-canopy ra, the profile, and the fixed 50. Integrating A over L = 0–2 gives crown photosynthesis, µmol per m² of ground per s:

| Tair | VPD | top ra for every layer (today) | doc profile | fixed 50 |
|---:|---:|---:|---:|---:|
| 15 | 1 | 29.5 | 31.0 (+5 %) | 29.8 |
| 25 | 1 | 28.8 | 28.7 (0 %) | 27.8 |
| 35 | 1 | 11.0 | 13.8 (+26 %) | 12.1 |
| 15 | 2 | 28.2 | 28.6 (+1 %) | 28.1 |
| 25 | 2 | 26.8 | 26.2 (−2 %) | 25.6 |
| 35 | 2 | 10.2 | 13.2 (+29 %) | 11.1 |

The profile's large gains at 35 °C are an artefact. At L = 1.5–2 the profile's ra is 206–385 s m⁻¹. Deep leaves there have little radiation to lose, but their transpiration is set by hydraulics, not by the boundary layer. The explicit balance then puts them 4–15 °C *below* air: 19.9 °C at Tair 35, which is under the wet-bulb temperature of 29–32 °C. A real leaf cannot do that, because a resistance that large would also throttle its vapour loss. Four of the 30 profile cells, and 8–9 of the 64 factorial cells at ra ≥ 110, sit below wet bulb. None does over the wind model's own range of 20–63 s m⁻¹. This is phylloptim #122 (the boundary layer resists heat but not vapour) showing up as a temperature.

**At stand level** (SCM, PM on, wind driver changed). U0 = 0.8 m s⁻¹ is ra = 50; U0 = 5 m s⁻¹ is ra = 20.

| Tair | one-year above-ground mass, kg m⁻² (U0 0.8 / 2 / 5; Fick) | offspring production (U0 0.8 / 2 / 5; Fick) | canopy-top Tleaf − Tair, °C (U0 0.8 / 2 / 5) |
|---:|---|---|---|
| 10 | 0.053 / 0.098 / 0.118; 0.102 | 2.2e-4 / 1.39 / 31.6; 100 | 10.0 / 6.7 / 4.4 |
| 15 | 0.026 / 0.083 / 0.142; 0.234 | 5.1e-7 / 5.4e-3 / 13.0; 311 | 9.6 / 6.4 / 4.2 |
| 20 | 0.0059 / 0.037 / 0.101; 0.339 | 8.6e-11 / 7.6e-6 / 0.30; 449 | 8.9 / 6.2 / 4.1 |
| 25 | 0 / 0.0055 / 0.031; 0.302 | 0 / 2.0e-10 / 2.8e-5; 388 | 7.5 / 6.1 / 4.1 |

**Outcome against the doc's 7.4 criteria.** Do not drop wind dependence: ra is a first-order input, and wind should be a real driver rather than a constant. Do not add the exponential profile: as specified, it needs a vapour boundary-layer conductance and a free-convection floor on ra first (phylloptim #122, #7).

## 7.5.3 — Net radiation

**What Rn is.** The radiation a leaf absorbs minus what it emits, in W m⁻²: the heat it must lose by evaporation or to the air. plant uses Rn = 2·PAR/4.57 − 40: shortwave taken as twice PAR, and a fixed longwave loss of 40 W m⁻². Rn and ra enter the leaf temperature as a product, so an error in either scales the warming.

The comparison is against the default Rn, with the same conventions as above:

| variant | change to Rn | cells with \|ΔA\| > 5 % (of 36) | median \|ΔA\| | median \|ΔTleaf\| (°C) | PAR 1000 / 25 °C / 2 kPa: Tleaf, A |
|---|---|---:|---:|---:|---|
| Rn × 0.7 | −30 % | 19 | 5.6 % | 1.6 | 27.9 °C, 14.6 |
| default | — | — | — | — | 31.0 °C, 11.5 |
| Rn × 1.3 | +30 % | 22 | 7.0 % | 3.1 | 34.2 °C, 8.1 |
| longwave at εa = 0.9 | −38 to −53 W m⁻², in place of −40 | 2 | 0.5 % | 0.1 | 30.9 °C, 11.6 |
| longwave at εa = 0.8 | −76 to −106 | 21 | 6.7 % | 1.2 | 29.8 °C, 12.8 |
| longwave at εa = 0.7 | −114 to −159 | 23 | 12.0 % | 2.9 | 28.6 °C, 13.9 |

The longwave rows replace the fixed −40 with the doc's own clear-sky form, 0.97·σ·Ta⁴·(εa − 1). Here εa is the sky's emissivity: about 0.7 for a clear dry sky, 0.9 for a humid or cloudy one. The range spans Tair 15–40 °C.

**Outcome against the doc's criteria.** The doc's rule is to simplify Rn if a ±30 % change moves outputs by less than 5 %. It moves A by more than 5 % in about half the cells, so Rn needs more care, not less. The longwave rule is to drop the term (and fix it at −40) if it varies by less than about 10 W m⁻² across realistic skies, or moves Tleaf by less than 0.5 °C. It varies by 76–106 W m⁻² across εa 0.7–0.9 and moves Tleaf by a median 1.2–2.9 °C. So keep the longwave term and model it: phylloptim #28, the temperature-dependent outgoing longwave this test was gating. The fixed −40 is the εa ≈ 0.9 sky, so under clear skies it overstates leaf heating by 1–3 °C (up to 5.8).

## 7.5.4 — Stand level

**What the quantities mean.** *One-year above-ground mass* is what a patch builds in its first year (kg m⁻²); it stands in for annual carbon gain. *Water uptake* is the water the stand drew from the soil in that year, against 1 unit of rainfall. *Offspring production* is TF24's fitness readout: offspring produced over a 5-year patch. *Canopy-top Tleaf − Tair* is the tallest cohort's leaf warming at the end of year 1.

| Tair (°C) | canopy-top Tleaf − Tair | one-year above-ground mass, PM vs Fick | one-year water uptake | 5-year above-ground mass | offspring, PM ÷ Fick (PM; Fick) | lma-0.0825 share of two-species offspring, Fick → PM |
|---:|---:|---:|---:|---:|---|---|
| 10 | +6.7 | −4 % | +34 % | −21 % | 0.014 (1.39; 100) | 0.91 → 0.91 |
| 15 | +6.4 | −65 % | −34 % | −33 % | 1.7e-5 (0.0054; 311) | 0.86 → 0.79 |
| 20 | +6.2 | −89 % | −70 % | −49 % | 1.7e-8 (7.6e-6; 449) | 0.84 → 0.73 |
| 25 | +6.1 | −98 % | −91 % | −68 % | 5e-13 (2.0e-10; 388) | 0.86 → 0.64 |
| 30 | +6.6 | −100 % (stand dies) | −100 % | −100 % | 0 (0; 112) | 0.96 → none |
| 35 | +7.9 | −100 % | −100 % | −100 % | 0 (0; 9.9e-6) | 1.00 → none |

**Why the stand collapses: the deficit, not the temperature.** Two Fick runs were matched to PM's canopy-top leaf. The first sets Fick's leaf temperature to PM's leaf temperature, keeping the air deficit of 1 kPa. The second also sets the deficit to the leaf-to-air deficit PM's leaf sees.

| Tair | PM: mass, offspring | Fick at PM's leaf temperature | Fick at PM's leaf temperature and leaf deficit | that deficit (kPa) |
|---:|---|---|---|---:|
| 10 | 0.098, 1.39 | 0.278, 372 | 0.097, 0.65 | 1.67 |
| 15 | 0.083, 0.0054 | 0.347, 455 | 0.081, 0.19 | 1.85 |
| 20 | 0.037, 7.6e-6 | 0.265, 335 | 0.037, 0.0022 | 2.07 |
| 25 | 0.0055, 2.0e-10 | 0.076, 54 | 0.0065, 1.2e-8 | 2.35 |

Warmth alone barely hurts: Fick at PM's leaf temperature builds 0.08–0.35 kg m⁻² and reproduces. Adding the leaf-to-air deficit reproduces PM's one-year growth to within 3 % at 10–20 °C (18 % at 25 °C). The mechanism shows at the canopy top. The PM leaf's photosynthesis is within 4–30 % of the Fick leaf's. But a leaf at 1.7–2.4× the deficit loses that much more water per unit of CO₂ it takes in, and TF24 charges a steeply rising hydraulic cost for the water. The PM leaf's stem water potential is 0.8–1.1 MPa more negative (in magnitude, 2.1 → 3.1 MPa at 15 °C), and the tallest cohort's net carbon gain falls by 10 % at 10 °C and by 71–99 % at 15–25 °C. TF24's default drivers were set up under Fick, where the leaf sees 1 kPa. PM raises the deficit its leaves actually face, and TF24 is very sensitive to that.

**Why offspring moves so much more than growth.** In these 5-year patches the tallest Fick plants reach 5.3–5.4 m, just past the 5 m height at which TF24 starts allocating to reproduction. PM plants reach 5.1 m at 10 °C and stay below 5 m from 15 °C up. A 21 % loss of 5-year mass at 10 °C therefore costs 99 % of the offspring. Read offspring here as "reproduces or not", not as a percentage.

**Competition.** The lower-lma species wins in every case where a stand survives. PM narrows its lead (from 86 % to 64 % of offspring at 25 °C) but does not reverse it anywhere. At 15 °C and above PM's absolute offspring are 10⁻² or less, so the competitive outcome under PM is moot.

**Outcome against the doc's criterion** (a change of more than 5 % in annual carbon gain under representative summer conditions): one-year growth moves by more than 5 % at every air temperature except 10 °C, and by 89–100 % at 20–35 °C.

## Caveats

- **Constant drivers.** TF24's SCM holds every driver constant: midday radiation (PPFD 1800) around the clock, 1 kPa, 2 m s⁻¹. "Annual" here means a year of constant midday conditions, not a seasonal or diurnal climate at a dry-sclerophyll site, as the doc envisaged. *Inferred*: a diurnal cycle would reduce PM's effect, since leaves heat least when radiation is low.
- **TF24 has not been recalibrated for PM.** The collapse says PM-on TF24 at today's defaults is not a viable stand model. It does not say PM is wrong.
- **The leaf runs use Step 5's tuned leaf** (root carbon 20 per unit leaf area, conductance 5 × 10⁻³), not a TF24 plant.
- **One boundary layer, for heat only** (phylloptim #122). There is no vapour boundary layer and no free convection (#7), so at large ra a transpiring leaf can fall below wet bulb. Over the wind model's own range of 20–63 s m⁻¹ none did.
- **The SCM runs use plant's solver as it is.** The stand-level leaves were not scanned for the missed-maximum problem. *Inferred*: at leaf level, at Tair ≤ 25 °C, the solver was within 0.05 µmol m⁻² s⁻¹ of the global maximum in every cell, so the 10–25 °C collapse is not an artefact of it. The misses are at Tair ≥ 35 °C and at the zero-flow end of the bracket.
- **Cost.** Each SCM scenario takes 1–7 s, except the 35 °C PM-on stand, which takes 10 minutes while it dies (61 minutes on a heavily loaded machine; the output was bit-identical).
