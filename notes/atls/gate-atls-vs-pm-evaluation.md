# ATLS verification gate — thermal damage feedback vs PM-only (#566)

The #566 evaluation gate, mirroring the #523 Step-5 factorial: does the
Lumry–Eyring leaf thermal-damage feedback (TF24t) **materially** change carbon
gain relative to PM-only (TF24 with the Penman–Monteith energy balance, no
damage), under representative Australian heatwave leaf temperatures — and does it
do so in the predicted place? Companion to `implementation-plan.md` (Phases 1–5
built the layer; this is the "keep only if material" measurement).

Regenerate: `Rscript notes/atls/gate_atls_vs_pm.R notes/atls/gate_atls_vs_pm`
(writes `gate_atls_vs_pm.rds` + `_A.csv`/`_B.csv`). The script is self-contained.

## Method

Four measurements, deliberately separating the *mechanism* from the two things
that otherwise confound it — the minimal-cut PM's leaf-overheating bias, and
TF24t's cost machinery:

- **(A) Clean mechanism** — a single TF24-default `Leaf` (Step-5 geometry) with
  the damage layer OFF vs ON, at a **prescribed** leaf temperature (Fick path,
  `Tleaf = leaf_temp`, no energy-balance overheating), swept `Tleaf ∈ 25…46 °C`.
  This isolates the pure `jmax·N` damage response.
- **(B) Realistic PM operating point** — the same leaf with `use_energy_balance`
  on, so N is evaluated at the midday operating-point `Tleaf` (the Phase-4
  wiring), swept over `Tair × PAR`. Honest about the minimal-cut PM.
- **(C) Whole-plant net production** (`net_mass_production_dt`, fixed height 8 m)
  across a prescribed-midday `Tair` gradient: TF24 (PM-only) vs **TF24t with all
  thermal costs zeroed** (isolates the damage feedback) vs full TF24t.
- **(D) Acclimation buffering** — at a prescribed `Tleaf = 42 °C` (within reach
  of `T_crit`'s acclimation ceiling `38 + 6 = 44 °C`), sweep the `A_crit` state
  and read N and A back.

All defaults are the ATLS/TF24t defaults (`T_crit_0 = 38 °C`, switch slope
`m = 1`, `k_d1_0 = k_r1_0 = 864 d⁻¹`, repair cut `t_rep_cut = 45 °C`).

## Results

### (A) The damage mechanism is real and large in the 36–44 °C band

| Tleaf (°C) |   N   | A off | A on | ΔA %  |
|-----------:|------:|------:|-----:|------:|
| 30         | 1.000 | 15.49 | 15.49|  −0.0 |
| 34         | 0.982 | 11.82 | 11.74|  −0.6 |
| 36         | 0.891 |  6.68 |  6.36|  −4.8 |
| 38         | 0.653 |  4.25 |  3.41| −19.9 |
| 40         | 0.500 |  2.68 |  1.68| −37.4 |
| 42         | 0.439 |  1.65 |  0.89| −46.3 |
| 44         | 0.375 |  0.99 |  0.40| (A<1) |
| 46         | 0.286 |  0.57 |  0.10| (A<1) |

(PAR = 1000; PAR = 2000 is within ~1 point.) N is a smooth turn-off centred just
above `T_crit = 38 °C`; assimilation loss reaches **46 %** by 42 °C. Below 34 °C
the layer is inert (N ≈ 1), so it does not perturb conditions where damage is
physiologically absent.

### (B) On the PM path the leaf overheats into the damage zone early

The minimal-cut PM heats the leaf **5–22 °C above air** (operating-point
`Tleaf`: 42 °C at `Tair = 32`, 48–60 °C at `Tair = 38`). So the damage zone is
reached at surprisingly *low* air temperatures — at `Tair = 32 °C, PAR = 1000`
the operating leaf sits at 42 °C and A drops **−52 %** vs PM-only. By `Tair ≥ 38`
the PM-only leaf has already shut down (A ≈ 0, ΔA undefined), so the *incremental*
damage there is small — the shutdown, not the damage, dominates the extreme
corner. This is the same leaf-overheating behaviour the #523 Step-5 gate flagged;
here it means the damage feedback bites well within the normal operating range.

### (C) Whole-plant: the effect is the damage feedback, not the costs

Net production (mg-equiv, height 8 m), relative change vs TF24 PM-only:

| Tair (°C) | net PM-only | Δdamage-only % | Δfull TF24t % |
|----------:|------------:|---------------:|--------------:|
| 25        |     +0.756  |         −1.4   |        −2.0   |
| 32        |     −3.501  |        −20.7   |       −21.5   |
| 38        |     −4.862  |         −6.8   |        −6.9   |
| 42        |     −5.112  |         −5.2   |        −5.6   |
| 46        |     −5.186  |         −2.1   |        −2.4   |

The **damage-only** arm (all cost coefficients zeroed) already accounts for
almost the entire gap; adding the full cost machinery moves it by <1 point. So
the whole-plant signal is the `jmax·N` feedback, not the ATLS respiration/
construction costs. The effect exceeds 5 % at `Tair` 32–42 °C. (Net production is
negative for `Tair ≥ 32`: a fixed 8 m plant under constant year-round heat
respires more than it fixes — the *relative* damage effect is the signal, not the
sign, which is an artefact of the constant-climate isolated-plant setup.)

### (D) T_crit acclimation buffers the damage

At a prescribed 42 °C leaf, raising the `A_crit` acclimation state recovers both
N and A:

| A_crit | N     | A     |
|-------:|------:|------:|
| 0.0    | 0.439 | 0.898 |
| 1.0    | 0.513 | 1.114 |
| 2.0    | 0.606 | 1.359 |
| 5.0    | 0.741 | 1.570 |

N climbs from 0.44 to 0.74 (A: +75 %) as `T_crit` acclimates from 38 toward its
44 °C ceiling — acclimation is a working, valuable strategy axis, exactly as
designed.

## Verdict — KEEP the ATLS damage layer

All sub-gates clear, and the effect is localised where the theory predicts:

- **Material leaf response** — ΔA up to −46 % in the 38–44 °C damage band, far
  above the 5 % bar. ✅
- **Material whole-plant response** — damage-only net production changes 5–21 %
  at `Tair` 32–42 °C. ✅
- **Inert where it should be** — N ≈ 1 and ΔA < 1 % below ~34 °C leaf
  temperature, so TF24t does not perturb non-heat conditions. ✅
- **Acclimation buffers it** — the `A_crit` axis measurably restores carbon,
  confirming the strategy space is meaningful. ✅

## Caveats / what this does NOT establish

- **Isolated leaf / fixed-size plant, tuned geometry.** Absolute magnitudes
  depend on the chosen conductance/root mass and on the constant-climate driver;
  the *pattern* (steep 36–44 °C damage band, inert below, acclimation recovery)
  is the robust result. A full SCM competitive-outcome run (does a tolerant/
  acclimating strategy *win* under an Australian heatwave climate vs a PM-only
  resident?) is the next evaluation and is not run here.
- **Entangled with the minimal-cut PM's leaf overheating.** Because PM heats the
  leaf 5–22 °C above air (`Rn = 2·PAR` with no NIR split, fixed `ra`/longwave,
  leaf-to-air VPD deferred — the #523 caveats), the *air* temperature at which
  damage engages is biased low. The quantitative gate therefore rests on the
  clean prescribed-`Tleaf` arm (A/D); the PM-path arm (B) shows the mechanism is
  reached in normal operation but should be re-measured once PM gains a
  leaf-to-air VPD feedback.
- **Cost coefficients are calibration targets**, not measured values (Phase 3).
  (C) shows the current defaults are second-order to the damage feedback, so the
  KEEP decision does not hinge on them.
- **Acclimation timescale debt (Phase 2).** The `A_crit` sweep in (D) is a static
  read at fixed states; the kinetics are labelled d⁻¹ but integrated on the
  yearly SCM clock, so the *rate* at which acclimation is reached in a run — and
  the induction-cost scale — awaits the timescale reconciliation.
