# Implementation Plan: Penman-Monteith Transpiration in `plant`

## 1. Background

The `plant` model currently computes transpiration using a simple Fick's law approximation:

```
E = gs · VPD
```

where `gs` is stomatal conductance (output of the stomata optimisation) and `VPD` is vapour pressure deficit of the ambient air. This is computationally cheap and sufficient when leaf temperature is close to air temperature, but it has two related limitations:

1. **Leaf temperature is assumed equal to air temperature.** In reality, the latent heat of transpiration cools the leaf, and high radiation loads can warm it. Under hot, high-radiation conditions typical of Australian summers, Tleaf can deviate 5–10°C from Tair.

2. **VPD at the leaf surface is wrong.** The driving gradient for transpiration is the difference between saturation vapour pressure at Tleaf and actual vapour pressure of the air. Using Tair to compute saturation vapour pressure introduces systematic error when Tleaf ≠ Tair.

Both errors propagate into the stomata optimisation: the water cost (E) is miscalculated, and the carbon benefit (A via Farquhar) is computed at the wrong temperature.

---

## 2. Options Considered

### Option 1: Keep Fick's law (status quo)
Simple and fast. Adequate when Tleaf ≈ Tair — i.e., cool conditions, low radiation, or where temperature-driven errors are secondary to other uncertainties. Not adequate for hot/dry Australian conditions.

### Option 2: Full iterative leaf energy balance
At each evaluation of the profit function: guess Tleaf → compute E → compute sensible heat H → update Tleaf → iterate until convergence. Physically exact, works at single-leaf scale. Used in MAESPA and similar detailed models. Adds iteration overhead inside the optimisation loop — a concern for `plant`'s performance-sensitive architecture.

### Option 3: Penman-Monteith (preferred)
Replace Fick's law with the Penman-Monteith equation. PM solves the coupled energy balance and vapour diffusion problem algebraically, eliminating the need for iterative solution by linearising the saturation vapour pressure curve around Tair. Gives physically consistent E and Tleaf simultaneously with no iteration. Adds modest new parameter requirements (Rn, ra) but no architectural change to the optimisation loop.

---

## 3. Preferred Path: Replace Fick's Law with Penman-Monteith

### 3.1 The Penman-Monteith Equation

The PM equation for latent heat flux at the leaf scale is:

```
λE = (Δ · Rn + ρ·cp · VPD / ra) / (Δ + γ · (1 + rs/ra))
```

where:

| Symbol | Meaning | Units |
|--------|---------|-------|
| λE | Latent heat flux (transpiration) | W m⁻² |
| Δ | Slope of saturation vapour pressure curve at Tair | Pa K⁻¹ |
| Rn | Net radiation at the leaf | W m⁻² |
| ρ·cp | Volumetric heat capacity of air | J m⁻³ K⁻¹ |
| VPD | Vapour pressure deficit (ambient air) | Pa |
| ra | Aerodynamic resistance (boundary layer) | s m⁻¹ |
| γ | Psychrometric constant | Pa K⁻¹ |
| rs | Stomatal resistance = 1/gs | s m⁻¹ |

Transpiration rate E (mol m⁻² s⁻¹) is recovered from λE by dividing by the latent heat of vaporisation λ and molar mass of water.

Leaf temperature follows directly from the sensible heat flux:

```
H = Rn - λE
Tleaf = Tair + H · ra / (ρ·cp)
```

### 3.2 Change to the Optimisation Loop

The only change required is to the E submodel called inside the profit function evaluation. Currently:

```
E = gs · VPD          # Fick's law
```

Replaced by:

```
rs = 1 / gs
λE = PM(Rn, VPD, ra, rs, Tair, P)   # Penman-Monteith
E = λE / λ                            # convert to mol m⁻² s⁻¹
Tleaf = Tair + (Rn - λE) · ra / (ρ·cp)
```

Tleaf is then passed to the Farquhar photosynthesis model in place of Tair. No other changes to the optimisation structure are required.

### 3.3 Net Radiation at the Leaf

Rn must be computed for each leaf position in the canopy. It has two components:

**Shortwave (absorbed solar):**
```
SW_abs = SW_inc · (1 - α_sw)
```

`plant` already tracks absorbed PAR. Since PAR comprises approximately 50% of total shortwave radiation, a first approximation is:

```
SW_abs ≈ 2 · PAR_abs
```

A more careful treatment would apply separate albedo values to PAR and near-infrared bands (leaves are more reflective in NIR: α_PAR ≈ 0.1, α_NIR ≈ 0.5).

**Longwave (thermal):**
```
LW_net = εa·σ·Ta⁴ - εl·σ·Tleaf⁴
```

This depends on Tleaf — which is what we are solving for. PM's linearisation handles this implicitly: it approximates LW_net ≈ f(Ta) using the slope Δ, avoiding the need to know Tleaf when computing Rn. For consistency, compute Rn using Tleaf = Tair (small error, typically < 5 W m⁻²) and rely on PM's linearisation for the remainder.

In practice:
```
LW_net ≈ ε·σ·(εa·Ta⁴ - Ta⁴) = ε·σ·Ta⁴·(εa - 1)
```

Since εa < 1 (typically 0.7–0.9 depending on cloud cover and humidity), LW_net is a net cooling term. A simple approximation is LW_net ≈ −30 to −50 W m⁻² under clear-sky conditions, which can be parameterised as a fixed offset or derived from Ta and humidity.

---

## 4. Wind Profile Submodel

Aerodynamic resistance `ra` is required by PM and varies with wind speed and canopy position. A simple but physically grounded approach is an exponential wind attenuation model.

### 4.1 Model

Above-canopy wind speed U₀ (at canopy top, height H) attenuates exponentially with cumulative LAI from the top:

```
U(L) = U₀ · exp(-α_w · L)
```

where L is cumulative LAI from the top of the canopy (0 at top, LAI_total at ground) and α_w is the wind extinction coefficient (typically 2.5–3.5 for closed canopies).

Aerodynamic resistance at each layer follows from leaf-level boundary layer theory:

```
ra(L) = C_ra · sqrt(d / U(L))
```

where d is characteristic leaf dimension (m) and C_ra is a constant (~200 s⁰·⁵ m⁻¹ in standard formulations).

This can be simplified to:

```
ra(L) = ra0 · exp(α_w · L / 2)
```

where ra0 is the aerodynamic resistance at the top of the canopy, computed from U₀ and d.

### 4.2 Integration with `plant`'s Canopy Structure

`plant` tracks vertical canopy structure via cumulative LAI. At each canopy layer (height z, cumulative LAI = L):

1. Compute U(L) from exponential profile
2. Compute ra(L) from leaf size and U(L)
3. Pass ra(L) to PM alongside layer-specific Rn and gs

This means ra varies continuously with canopy depth — deeper, more shaded leaves have higher ra (lower wind), which reduces their aerodynamic conductance and makes their transpiration more radiation-limited relative to demand-limited.

---

## 5. Implementation Steps

### Step 1: Implement PM as a standalone function
Write and test `transpiration_penman_monteith(Rn, VPD, ra, gs, Tair, P)` returning λE, E, and Tleaf. Unit test against known solutions.

### Step 2: Implement Rn submodel
Write `net_radiation_leaf(PAR_abs, Tair, ea, epsilon_a)` returning Rn. Initially use the PAR × 2 shortwave approximation and a fixed or humidity-derived LW_net term.

### Step 3: Implement wind profile submodel
Write `wind_speed_canopy(U0, L, alpha_w)` and `aerodynamic_resistance(U, d)`. Validate that ra values are physiologically reasonable (~20–200 s m⁻¹).

### Step 4: Wire into optimisation loop
Replace `E = gs * VPD` with the PM call. Pass Tleaf to the Farquhar model. Ensure all new inputs (U0, d, PAR_abs, Rn components) are threaded through from the environment object.

### Step 5: Sensitivity analysis
Run factorial comparisons across a range of conditions (Tair: 15–40°C, VPD: 0.5–4 kPa, PAR: 200–2000 µmol m⁻² s⁻¹) comparing Fick's law vs PM outputs for E, Tleaf, and A. Identify the conditions under which the correction is most material.

### Step 6: Validate against leaf-level measurements
Where possible, compare PM-predicted Tleaf against thermal infrared measurements from field campaigns in the AusTraits or related datasets.

---

## 6. Parameters

### Fixed physical constants (hardcode)

| Parameter | Symbol | Value | Units | Notes |
|-----------|--------|-------|-------|-------|
| Latent heat of vaporisation | λ | 2.45 × 10⁶ | J kg⁻¹ | Weak function of T; fix at 25°C |
| Psychrometric constant | γ | 66.5 | Pa K⁻¹ | At sea level; scale with P |
| Stefan-Boltzmann constant | σ | 5.67 × 10⁻⁸ | W m⁻² K⁻⁴ | Physical constant |
| Leaf emissivity | εl | 0.97 | — | Near-universal for leaves |
| Volumetric heat capacity of air | ρ·cp | 1200 | J m⁻³ K⁻¹ | Weak function of T and P |

### Derivable from environment (compute at runtime)

| Parameter | Symbol | How to compute |
|-----------|--------|----------------|
| Slope of sat. vap. pressure curve | Δ | `Δ = 4098 · es(Tair) / (Tair + 237.3)²` (Tetens formula) |
| Saturation vapour pressure | es | `es = 0.6108 · exp(17.27·T / (T + 237.3))` kPa |
| Atmospheric emissivity | εa | `εa = 0.642 · (ea/Ta)^(1/7)` (Brutsaert 1975); or fix at 0.80 under overcast |

### New model parameters (to be set in species/environment object)

| Parameter | Symbol | Suggested default | Units | Sensitivity |
|-----------|--------|-------------------|-------|-------------|
| Leaf shortwave albedo | α_sw | 0.45 | — | Low; PAR→NIR ratio matters more |
| Leaf dimension (for ra) | d | 0.05 | m | Moderate; use species mean or guild value |
| Wind extinction coefficient | α_w | 2.5 | — | Moderate; well-constrained for closed canopies |
| Above-canopy wind speed | U₀ | 2.0 | m s⁻¹ | High; drive from environment if possible |
| Atmospheric emissivity (clear sky) | εa | 0.80 | — | Low for daytime Rn; larger effect on nocturnal budget |

### Parameters already in `plant` (no change needed)

- gs (output of stomata optimisation → rs = 1/gs)
- PAR absorbed (per canopy layer)
- Tair, VPD, atmospheric pressure P (environment drivers)
- Cumulative LAI per canopy layer (canopy structure)

---

## 7. Critical Assessment: Is the Added Complexity Justified?

### 7.1 The Overparameterisation Problem

A well-founded critique of land surface model development — articulated forcefully in the PLUMBER benchmarking work (Abramowitz et al.) — is that LSMs are systematically overparameterised relative to the observations available to constrain them. Each process addition is locally justified, but the aggregate model accumulates degrees of freedom that cannot be identified from flux tower or eddy covariance data. The consequence is equifinality: many parameter combinations fit available data equally well, model skill does not improve with complexity, and the appearance of physical realism masks unidentifiable assumptions.

This critique deserves to be taken seriously here. The PM addition introduces new parameters (U₀, d, α_w) that vary substantially in the real world but will in practice be fixed at assumed values — because the data to constrain them spatially and temporally does not exist at the scale `plant` operates. There is a real risk of adding the appearance of process fidelity without improving predictive skill or biological insight.

### 7.2 Where the Critique Applies to This Addition

The wind profile submodel is the most exposed component. U₀ varies by orders of magnitude across sites and through time; α_w depends on canopy architecture in ways that are hard to measure; and ra is sensitive to both. In practice these will be fixed or given weak priors, meaning the wind-driven variation in ra is an assumption, not a constraint. If the key model outputs are insensitive to ra across realistic ranges (20–200 s m⁻¹), then the wind submodel adds complexity without value and should be dropped — replacing it with a single fixed ra.

Atmospheric emissivity εa and the longwave Rn term are similarly uncertain. Cloud cover, humidity, and aerosol loading all affect εa in ways that `plant`'s environment forcing typically does not resolve. A fixed LW_net offset is arguably more honest than a formula that implies precision that isn't there.

### 7.3 Where the Addition Is Genuinely Justified

`plant` is not an LSM in the PLUMBER sense. Its goal is mechanistic understanding of plant growth and competition, not predictive skill at flux tower sites. This shifts the calculus in two ways.

First, the stomata optimisation already embeds a strong mechanistic prior — gs is derived, not fitted. PM is the physically consistent extension of that philosophy to the water flux calculation. Using Fick's law with Tleaf = Tair is also a modelling choice; it just makes its assumptions implicitly rather than explicitly.

Second, the Tleaf → Farquhar feedback is a *structural* correction, not parameter elaboration. Getting Tleaf wrong introduces directional bias into A: hotter leaves have lower Vcmax efficiency and higher respiration, so underestimating Tleaf under high radiation systematically overestimates carbon gain. This is not noise — it is a systematic error that will be largest in exactly the conditions (hot, dry, high radiation) most relevant to Australian vegetation modelling. The correction is well-motivated independent of whether the wind profile submodel is retained.

### 7.4 Decision Framework: What Justifies Keeping Each Component

The implementation should be treated as a staged hypothesis test, not a commitment. Each component is kept only if it demonstrably changes outputs in scientifically material ways.

**The PM core (replacing Fick's law with PM + Tleaf feedback to Farquhar)**

Keep if: the sensitivity analysis (Step 5) shows that Tleaf deviates meaningfully from Tair under target conditions (Tair > 30°C, PAR > 1000 µmol m⁻² s⁻¹, VPD > 2 kPa), and that this deviation materially changes predicted A, gs, or growth rates. A threshold of >5% change in annual carbon gain under representative Australian summer conditions is a reasonable bar.

Drop if: Tleaf − Tair < 2°C across realistic conditions, or the downstream effect on A is negligible. In that case Fick's law is adequate and simpler.

**The Rn submodel**

Keep the PAR × 2 shortwave approximation if: the sensitivity of PM outputs to Rn uncertainty is large enough that a factor-of-two error in the NIR component would matter. If PM outputs are dominated by the VPD/ra term rather than the Rn term (common when VPD is high), Rn precision matters little and a fixed or highly simplified Rn is adequate.

Drop the longwave component if: LW_net varies less than ~10 W m⁻² across realistic εa values, or if its effect on Tleaf is < 0.5°C. In that case fix LW_net at −40 W m⁻² and acknowledge the assumption.

**The wind profile submodel**

Keep if: ra varies enough with canopy depth to materially change the within-canopy distribution of E and Tleaf, and if this distribution matters for model outputs (e.g. canopy carbon gain, competitive dynamics between canopy layers). This requires comparing outputs with fixed ra vs. depth-varying ra.

Drop (replace with fixed ra) if: the outputs of interest are insensitive to within-canopy ra variation, or if U₀ uncertainty dominates ra uncertainty to the point where the profile model adds nothing over a single fixed value. A fixed ra of 50 s m⁻¹ is a reasonable fallback representing moderate wind conditions at mid-canopy.

### 7.5 Recommended Evaluation Protocol

Before committing to any component, run the following comparisons and report effect sizes explicitly:

1. **Fick vs PM** across factorial combinations of Tair (15, 25, 35, 40°C) × VPD (0.5, 1, 2, 3 kPa) × PAR (200, 500, 1000, 2000 µmol m⁻² s⁻¹). Report ΔTleaf, ΔE, ΔA, and Δgs. If the maximum ΔA across this space is < 5%, do not proceed.

2. **Fixed ra vs wind profile** across the same factorial, with ra fixed at 50 s m⁻¹ vs. the exponential profile. Report the range of ra produced by the profile model and its effect on E and Tleaf.

3. **Rn sensitivity**: vary Rn ± 30% (representing uncertainty in PAR→SW conversion and LW_net) and report sensitivity of PM outputs. If outputs change < 5%, simplify the Rn submodel.

4. **Annual integration**: run a full annual simulation under representative Australian conditions (e.g. a dry sclerophyll site) with Fick vs PM and report effect on annual carbon gain, water use, and competitive outcomes. This is the ultimate test of whether the addition matters for the science.

Components that pass their test are retained. Components that do not are dropped, with the null result documented explicitly. This prevents complexity accumulation and provides a defensible audit trail.

---

## 8. Priority and Scope

A minimal first implementation requires only Steps 1–4 and can use fixed values for d, α_w, and U₀. This delivers physically consistent E and Tleaf with three new parameters (d, α_w, U₀) that are well-constrained and weakly sensitive. The full treatment with sensitivity analysis and field validation (Steps 5–6) is recommended before using PM-derived outputs in publications, and the evaluation protocol in Section 7.5 must be completed before any component is treated as an improvement over the status quo.

The improvement is expected to be most material for:
- High-temperature, high-radiation conditions (Australian summer)
- Deep canopy layers where ra is large and the radiation-driven term dominates
- Species with large leaves (high d → high ra → stronger radiation coupling)
