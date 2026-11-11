# plant#275: can TF24 reproduce Ross et al. 2008's Callitris dynamics?

Run on plant develop@0b1a4b7a (TF24@v15), fieldtrial main@0c1179e, phylloptim e2c9460, on 2026-11-11. The draft set `callitris_central_west_nsw.yml` lives here until it moves to fieldtrial. WorldClim 2.1 at Ross's study-area centre (148°E, 33°S) gives 572 mm yr⁻¹, inside Ross's 470–595 mm. Traits are AusTraits geometric means for *C. glaucophylla* (austraits-means 7.0.0.9000): LMA 258 g m⁻², wood density 631 kg m⁻³, seed mass 9.8 mg, Huber value 4.7e-4. Height: `hmat` 12 m with `a_f2` 4.4, so reproductive allocation reaches 0.95 at 20 m (mid Ross's 16–23 m site maxima). This is a calibration proxy: AusTraits heights are maxima, and the maximum comes from `hmat` plus the allocation curve. Zotero, mulga-library and austraits-means hold no Callitris P50, so hydraulics are plant defaults, with an uncited "tolerant" variant (P50s and `psi_crit` ×2) where the default fails. Phenomena come from Ross-2008 (doi:10.1016/j.ecolmodel.2007.08.012) and the review it cites, Thompson-2005a (doi:10.1071/bt04115).

| # | Phenomenon (source) | TF24 mechanism | Captured? | Evidence |
|---|---|---|---|---|
| P0 | Long-lived adults; drought deaths rare (Thompson-2005a; Ross Fig. 6b, 19.1 of 20 survive 50 y) | `d_I` + storage-driven mortality; hydraulics | **No at defaults**; yes with parameter changes | Default hydraulics: a lone 20 ha⁻¹ stand falls to 0.02 by 60 y, while soil stays at ψ ≈ −0.02 MPa. Leaf ψ sits at stem P50 (2.9 MPa), so this is a height/pathway limit, not water supply. Tolerant variant: 9.6 of 20 survive 60 y (the `d_I` = 1 % yr⁻¹ floor). Adding `d_I` 0.001 and `a_dG1` ×0.1: 20.5 of 22 over 50 y, 16.5 m (Ross 19.1, 19.4 m). |
| P1 | Episodic recruitment in wet years (Ross §4.3, 1950/1955; Thompson-2005a waves) | Establishment = f(seedling carbon gain at germination), reading soil water; rain drivers or pulses | **Weakly** | Spline climatology: the open-stand soil never dries, and establishment is 0.997 every year, wet or not. Poisson storm rain: establishment ranges 0.25–1.0 between dry and wet quarters (0.03–0.8 with `a_d0` ×10), but recruits arrive every year; wet-year cohorts are at most ~2× their neighbours. |
| P2 | Canopy suppresses recruitment (Ross Fig. 4, none above 14 m² ha⁻¹; Fig. 7b, none in stands ≥ 5,000 ha⁻¹) | One-way light competition; water at establishment | **No** | Ground light never drops below ~0.6, because every stand self-thins to ~1,200 ha⁻¹ by 20 y. Recruits by 25 y are 30–41k ha⁻¹ under no canopy, 1,000 ha⁻¹ and 30,000 ha⁻¹ starts alike. |
| P3 | "Lock-up": dense regrowth stagnates, no self-thinning (Ross §3.2, 10k–100k ha⁻¹ at 3–5 m; Thompson-2005a) | Light + water competition; storage mortality | **No (structural)** | 30,000 ha⁻¹ → 5,800 (10 y) → 1,240 (20 y) → ~300 (50–60 y) at 14 m. 1k, 10k and 30k starts converge, with losses in summer die-offs as the soil dries to −5 to −9.5 MPa. Storage mortality ×0.1 only slows it (626 ha⁻¹ at 50 y, 12 m). Suppressed trees starve; none persists with near-zero growth. |
| P4 | Regrowth suppresses overstorey growth (Ross Fig. 6: 14.7 vs 19.4 m; 16.7 vs 19.1 trees) | Water only, since shading is one-way | **Direction yes, form wrong** | With recruits the overstorey is wiped out (312 → 0.5 ha⁻¹; 20.5 → 0) and its height is unchanged. Ross has height suppressed and survival barely touched. |
| P5 | Thinning: release, then re-fill within 5–10 y when ≤ 1,000 ha⁻¹; retained canopy blocks recruits (Ross §4.4, Table 2) | `harvest(fraction, size band)` | **Mechanics yes, response partly** | Thinning at 4 m raises height 20 y later only from 11.0 to 12.1 m (Ross ~5 → 10 m). Every treatment re-fills, the unthinned one included (25k vs 48k recruits by +5 y), so the contrast that matters is missing (P2). |
| P6 | Site moisture sets growth and maximum height (Ross §3.4, Fig. 3) | Rain driver, soil | **No (emergently)** | Rain ×0.8 to ×1.25 gives identical heights (15.2–15.3 m at 60 y) and survival of 183–275 ha⁻¹. Ross also imposes it through site maximum height, so in TF24 it is a per-site `hmat`/`a_f2`, i.e. set-level. |
| P7 | Fire kills small trees (< 1–2 m) and crowns do not recover (Thompson-2005a; AusTraits `fire_killed` ×8) | `harvest(size_max)` within a stand; patch disturbance for stand-replacing fire | **Yes** | Killing 95 % of trees under 1.5 m at 30 y removed 40 % of recruits and left the overstorey intact; an all-size fire removes 95 % of everything. Limits: one kill fraction per size band; post-fire seed depression (up to 7 y) needs a birth-rate driver, because seed rain is external to the patch. |
| P8 | Grazing suppresses recruits until they escape at ~1 m; rabbits halted regeneration until myxomatosis (Ross §3.2–3.3; Thompson-2005a) | Repeated `harvest(size_max = 1)` | **Yes, qualitatively** | Removing 50 % of recruits under 1 m every quarter for 10 y: 3.9k vs 34.8k recruits ha⁻¹ (−89 %), with the overstorey kept at 462 vs 89 ha⁻¹. Release gives a wave (51k vs 35k by +5 y). Browsing as 40 discrete events is clumsy; TF24 seedlings pass 1 m in 1–2 y, faster than Callitris. |

## What it would take

**New model work** (plant): a mortality response to suppression that stops growth rather than killing (P3, which also drives P2 and P4); a recruitment window that integrates seedling survival over the first year or two of soil water, or a moisture trigger for germination (P1); possibly seedlings competing for water with established trees, or deeper canopy closure (P2). Events would be cleaner with a continuous small-size mortality rate for browsing and size-dependent fire kill, and with seed supply coupled to the patch's own trees for canopy seedbank and post-fire seed depression (P7, P8).

**Set work** (fieldtrial): a Callitris vulnerability curve (P50s), adult background mortality, the height curve (`hmat`, `a_f2`), sandy red-earth soil, and a real rain record (AGCD for Forbes 1948–1998, as daily pulses; the archive is not on this machine).

## Recommendation

Keep #275 in plant as model work. The paper's defining dynamics, lock-up and canopy-limited, wet-year recruitment, are not reachable by reparameterising TF24. The disturbance side (fire, grazing, thinning) already works through events. Spin the set off to fieldtrial as a sub-issue, `callitris_central_west_nsw`. It is a prerequisite for testing the model work, not a solution, and its open values (P50, mortality, height curve, AGCD rain) are fieldtrial's to source.

## Files

Run each script from this folder with plant, fieldtrial and their pinned dependencies installed: `Rscript expA.R`. The first run downloads WorldClim's global grids (~200 MB) into fieldtrial's cache. Each script takes a minute or two; the `.rds` results they write are not committed.

| File | What it does |
|---|---|
| `callitris_central_west_nsw.yml` | Draft fieldtrial set (schema 1): site, AusTraits traits, the `hmat` proxy, WorldClim drivers |
| `helpers.R` | Builds runs from the set: identical-trait species as separate cohorts, seed pulses as variable birth rates on a regular node schedule, rain scaling by year, stand and recruit summaries. Defines the `HEIGHT` and `TOLERANT` parameter overrides |
| `diag.R`, `diag2.R`, `soil.R` | Why an open stand dies at default hydraulics: leaf water potential, temperature, Huber value, soil water by layer |
| `expA.R` | P0, P3, P6: single cohorts at 20–30,000 ha⁻¹ for 60 y, both hydraulic variants, rain ×0.8 and ×1.25 |
| `expB.R` | P1, P2, P4: overstorey sown at 0, recruits from year 20, wet years 22 and 27 |
| `expC.R` | P5, P7, P8: thinning at 4 m, fire as size-banded harvest, quarterly grazing then release |
| `expD.R` | P0, P3, P4 with lower mortality; Ross §3.6 historical scenario |
| `expE.R`, `expF.R` | P1, P2 with Poisson storm rain; `expF.R` also multiplies `a_d0` by 10 |

Two traps found on the way. `run_scm(collect = TRUE)` records history only at node introductions, so a seed pulse made by cutting the node schedule short leaves gaps that `integrate_over_size_distribution()` fills with phantom stems; `helpers.R` keeps one regular schedule and carries the pulse in the birth rate. And `make_initial_state()` starts TF24 trees with zero storage, i.e. maximal storage mortality, so stands here are grown from seed instead.

Caveats: one draft parameterisation; tolerant hydraulics, lower mortality and `a_d0` are uncited sensitivities; climatology or synthetic rain, not a record; SCM densities are deterministic where Ross's are stochastic; seed rain is external to the patch.
