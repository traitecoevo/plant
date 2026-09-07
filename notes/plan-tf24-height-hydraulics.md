# TF24: height-dependent stem hydraulic resistance

## Implementation plan

Status: draft for discussion. Epic [#615](https://github.com/traitecoevo/plant/issues/615), branch `feature/tf24-height-hydraulics`, based on `develop`. Sub-epic of [#424](https://github.com/traitecoevo/plant/issues/424). Australian spelling. All structural claims flagged for demonstration over parameter sweeps rather than fit.

---

## 1. Problem statement

TF24 currently imposes resistance increasing linearly with height — [`src/tf24_strategy.cpp:385`](../src/tf24_strategy.cpp#L385), `leaf_specific_conductance_max = pars.K_s * pars.theta / (height * eta_c)`. Over a seedling-to-mulga range (0.3 → 8 m) this gives a 27-fold resistance increase, which greatly limits productivity across the size range where most demographic action occurs.

An earlier draft of this plan rationalised the linear form as the correct asymptote for the tallest trees, reached once conduit widening saturates. **That framing is retracted.** It required a conduit-diameter cap whose mechanism is unknown (see §2.1), and it was appealing largely because it flattered the existing implementation. Absent a mechanistic plateau, linear scaling is a convenient default, not an asymptote.

The replacement is not a new empirical fit. It is a spatially resolved anatomical model, integrated to a closed form, in which the height exponent is *derived* from two independently measurable quantities (the widening exponent and the Huber-value profile) rather than being a free parameter.

### Naming — `theta` retained, definition narrowed to the tip

TF24/FF16 uses `theta` for leaf area per unit sapwood area, currently constant across the whole plant. **Keep the name; narrow the definition to the terminal segment.**

The justification is structural rather than convenient: **the tip is the only position whose θ is height-invariant.** θ at the base necessarily changes with H — that *is* the compensation mechanism — so basal θ is a state variable, not a trait. A trait must live at the tip, and `theta` is therefore the correct name for it rather than a compromise. Symbol `theta_tip` / $\theta_{tip}$ is used in this document for clarity where position matters, but the code parameter stays `theta`.

Only the cost term is renamed: use `hcost` / 𝒞 for normalised hydraulic cost, avoiding the conventional Sperry θ. That is the newer and less embedded symbol.

**Correction — there is nothing to rename in code.** Neither plant nor phylloptim ever spells the hydraulic cost `theta`; it is `hydraulic_cost_`, `hydraulic_cost_TF()`, `hydraulic_cost_Sperry()` throughout, and `hydraulic_cost_` is an R-visible `Leaf` field. Introducing `hcost` would be a *new* name for a quantity that already has a clear one. Keep 𝒞 as notation in this document and leave the code alone.

The live `theta` collision is a different one, and this document missed it: phylloptim's `Leaf$theta_`, `theta_w_` and `theta_fc_` are **soil water content, wilting point and field capacity** (m³ m⁻³) for the Medlyn soil-moisture stress term, and all three are exposed to R. So `s$pars$theta` (m² m⁻², leaf per sapwood) and `leaf$theta_` (m³ m⁻³, soil water) are reachable in one session under the same name. Phase 2a's `theta_c` joins that stem. Dimensionlessness is the only thing distinguishing it.

#### Silent-breakage hazard

At β = 0.4 with $H_a$ = 20 m and `L_min` = 0.02 m:

$$\theta_{tip} = \theta_{old}\times 1000^{0.4} \approx \mathbf{15.8\times}$$

Same name, same units (m² m⁻²), a number ~16× different. Nothing in the type system catches this, and stale parameter files, saved runs, published tables and students' scripts would all be silently wrong. Three safeguards, in order of reliability:

1. **Bump the parameter-file schema version.** Files without the new tag are *rejected*, not reinterpreted. This is the only safeguard that actually catches the error, since units are unchanged. ⚠️ **This machinery does not exist and must be built.** plant has `scientific_version` / `model_id()`, a *model* version consumed by logpile's cache — not a file-format gate. RcppR6 pars lists carry no version tag, there is no `Parameters` serialiser, and nothing anywhere can reject an old parameter file. I11 is therefore a build task at Phase 2b, not a bump.
2. **Range guard** — an old-style value arrives ~16× too small; warn on implausible tip Huber values.
3. **Conversion helper** in-repo, recording `H_anchor` alongside every converted value.

Timing is favourable: in Phase 2a, $a_\theta = 0$ and `theta` has no profile, so its meaning is unchanged. **The narrowing lands at 2b only**, together with the schema bump.

#### Derived output to report

Expose **`theta_base(H)`** as a derived model output. This is the quantity comparable to McDowell et al.'s Douglas-fir series, so P2 becomes a clean prediction check without `theta` doing double duty as both trait and observable.

---

## 2. Within-plant anatomical profiles

Coordinate: **L** = distance from apex along the flow path (m). **H** = plant height. **L_min** = length of the terminal segment (m), a lower cutoff.

### 2.1 Conduit diameter — widening

$$D(L) = D_{tip}\left(\frac{L}{L_{min}}\right)^{b}, \qquad b \approx 0.2$$

**No upper cap on D.** A $D_{max}$ term was considered and rejected. There is a real empirical pattern behind it — the plateau in widening reported for 86–105 m *Sequoia*, *Sequoiadendron* and *E. regnans* — but the mechanism is unknown, and none of the candidate justifications applies here: freezing embolism is irrelevant to mulga or lowland eucalypt; conduit-level safety–efficiency coupling is empirically weak; and a packing or lumen-fraction ceiling is already represented through $F$ rather than needing a separate diameter cap.

Including it would put an unexplained empirical regularity into a parameter that then does load-bearing work bounding ESS height — the failure mode §6.4 and I4 exist to catch. It also introduces a derivative discontinuity in $R_L(H)$, which in an adaptive dynamics setting generates singular strategies at the kink more or less mechanically.

The quantity is retained as a **diagnostic ceiling** instead (§4.2), not as a model constraint.

*Caveat to log:* tracheid diameters are far more tightly bounded than vessel diameters. If TF24 is ever applied to conifers, a plateau may bite at much lower heights. The §4.2 diagnostic catches this case automatically.

### 2.2 Conduit density — the packing limit

Conduit number is **not** a free profile. Under a conserved lumen fraction $F = n_A \pi D^2/4$:

$$n_A(L) = n_{A,tip}\left(\frac{D(L)}{D_{tip}}\right)^{-2} \propto L^{-2b}$$

Widening is paid for by proportionally fewer conduits. This is the single most consequential correction relative to a naive single-tube treatment.

### 2.3 Sapwood-specific conductivity — D² not D⁴

$$k_s(L) \propto n_A(L)\,\langle D^4\rangle \propto F\,D(L)^2$$

$$\boxed{k_s(L) = k_{s,tip}\left(\frac{L}{L_{min}}\right)^{2b}}$$

**The exponent on conductivity is 2b, not 4b.** The $D^4$ of Hagen–Poiseuille survives only along a single continuous conduit; at the tissue level the packing limit halves the effective exponent.

### 2.4 Sapwood area and the Huber profile

$$A_s(L) = \frac{A_L(L)}{\theta(L)}$$

where $A_L(L)$ is leaf area distal to $L$. The key structural change: **θ becomes a within-plant profile, not a constant.**

$$\theta(L) = \theta_{tip}\left(\frac{L}{L_{min}}\right)^{-a_\theta}, \qquad a_\theta = 2b \text{ (default)}$$

Rationale: under one-conduit-per-leaf sectoriality, sapwood area per unit distal leaf area must grow as $D^2$ to accommodate wider conduits. Setting $a_\theta = 2b$ makes this a derived consequence with **zero new fitted parameters**, and yields Huber value ($1/\theta$) increasing basipetally as $L^{2b}$.

Keep $a_\theta$ as an explicit switch/sweep parameter rather than hard-coding $2b$ — see §7.

---

## 3. Trait set

| Symbol | **Code name** | Meaning | Units | Default | Evolvable? | Measurement |
|---|---|---|---|---|---|---|
| `ks_tip` | **`K_s`** (unchanged) | sapwood-specific conductivity, terminal segment | kg m⁻¹ s⁻¹ MPa⁻¹ | reparameterise from current | **yes** | twig segment conductivity |
| `theta` | **`theta`** (unchanged) | leaf area per sapwood area, **terminal segment** (narrowed definition) | m² m⁻² | reparameterise from current, ~16× larger | **yes** | terminal shoot dissection |
| `b` | **`D_c`** | widening exponent | – | 0.20 | no (conserved) | tip-to-base conduit series |
| `a_theta` | **`theta_c`** | Huber profile exponent | – | 2·`D_c` = 0.40 | no (derived) | derived; check vs data |
| `L_min` | **`L_tip`** | terminal segment length | m | 0.02 | no | **must match `K_s` protocol** |
| `D_tip` | — (not a parameter until §4.2) | tip conduit diameter | µm | species | no | anatomy; needed only for the §4.2 diagnostic |
| ρg | — | 0.00979 | MPa m⁻¹ | constant | no | physical constant |

**The code names are the contract.** `b` was unusable: `TF24_Pars` already has `b`, the Weibull vulnerability *scale*, derived in-struct from `p_50` and again in `make_TF24_hyperpar`. The `_c` suffix means "exponent" in this codebase (`c` and `root_c` are both curve shape exponents). `K_s` is not renamed to `ks_tip`: it is an R-facing *trait* read by the hyperpar, so renaming would break every `trait_matrix(..., "K_s")` call for a suffix a comment conveys — its meaning narrows at Phase 2a, its name does not move.

Two costs of `theta_c`, accepted rather than fixed: phylloptim spells soil water content `theta_`, `theta_w_`, `theta_fc_`, all R-visible, so `s$pars$theta_c` (dimensionless) and `leaf$theta_` (m³ m⁻³) coexist in one session; and `D_c` names a diameter profile the model never evaluates — `D_c` reaches the model only through β.

### 3.1 Identifiability hazard — flag prominently

For $H \gg L_{min}$:

$$R_L \propto \frac{\theta_{tip}}{k_{s,tip}}\,L_{min}^{\,2b+a_\theta}\,H^{\,1-2b-a_\theta}$$

The combination $\theta_{tip}\,L_{min}^{2b+a_\theta} / k_{s,tip}$ is the **only identifiable grouping**. `L_min` is therefore not an innocuous numerical cutoff — it enters with elasticity $2b + a_\theta \approx 0.6$ and trades off exactly against both `ks_tip` and `theta`.

**Protection rule:** `ks_tip`, `theta` and `L_min` must all be sourced from a **single terminal-segment definition** — the segment over which tip conductivity and tip leaf:sapwood area were actually measured. None may be calibrated independently of the others.

### 3.2 Reparameterisation from current TF24 values

Existing whole-plant `theta` corresponds to the *base* of a plant at some implicit anchor height:

$$\theta_{tip} = \theta_{current}\left(\frac{H_{anchor}}{L_{min}}\right)^{a_\theta}$$

**The anchor height is a real modelling decision with demographic consequences, not a bookkeeping detail.** Anchoring at seedling stature makes all taller plants cheaper; anchoring at mature stature makes seedlings more expensive than at present, shifting establishment and hence emergent size distributions. Anchor where the measurements are (mature), and treat `H_anchor` as a swept quantity.

See **§4.4** for the exact back-compatible parameter combination and the consistency check it yields.

---

## 4. Total resistance — closed form

$$R_L(H) = A_L\int_{L_{min}}^{H}\frac{dL}{k_s(L)\,A_s(L)}$$

For a trunk-dominated plant (leaf area concentrated in a crown thin relative to H), with $\beta \equiv 2b + a_\theta$:

$$\boxed{\;R_L(H) = \frac{\theta_{tip}\,L_{min}}{k_{s,tip}}\cdot\frac{(H/L_{min})^{\,1-\beta}-1}{1-\beta}\;}$$

Use the logarithmic limit $R_L = (\theta_{tip}L_{min}/k_{s,tip})\ln(H/L_{min})$ when $\beta \to 1$.

Note the sapwood-area *profile* drops out of the exponent entirely — it appears only in the prefactor, because $A_s$ enters both the integrand and the leaf-area normalisation.

### 4.1 Regimes of the height exponent

| $a_\theta$ | β | Height exponent | Interpretation |
|---|---|---|---|
| 0 (constant θ — current FF16) | 0.4 | **0.60** | widening only |
| 2b (derived, recommended) | 0.8 | **0.20** | widening + Huber compensation |
| > 1−2b | >1 | saturating | resistance approaches finite asymptote |

Over 0.3 → 8 m: linear = 27×, $H^{0.6}$ = 7.2×, $H^{0.2}$ = 1.9×.

**Correction — those are the asymptotic ratios, and they understate the closed form.** Evaluated on the actual expression with `L_tip` = 0.02 m over the same 0.3 → 8 m:

| β | asymptotic (above) | closed form |
|---|---|---|
| 0 | 27× | **28.8×** |
| 0.4 | 7.2× | **8.8×** |
| 0.8 | 1.9× | **3.3×** |

The gap is the $-L_{tip}^{1-\beta}$ term, which bites hardest at seedling stature — exactly where P6 lives. Use the closed-form column when arguing P5 or P7; the asymptotic one flatters the compensation by nearly 2× at β = 0.8.

**Consistency constraint to state explicitly in the paper:** you cannot independently specify `b`, the θ profile, and the resistance exponent. Pick two, derive the third. Constant θ combined with an assumed $H^{0.2}$ would double-count the compensation.

### 4.2 Diagnostic ceiling — falsification, not constraint

Pure $H^{1-\beta}$ applies across the whole height range, with no saturation regime. Justification: any plateau would only bite around 60–80 m (backing out the prefactor at b = 0.2 from Petit et al.'s measured basal hydraulic diameters for *E. regnans* of 220–250 µm), mulga never approaches it, and Petit et al.'s over-compensation result suggests widening in *regnans* is if anything stronger than needed rather than saturating early.

The maximum conduit diameter is instead used as a **post hoc check**:

$$D_{implied} = D_{tip}\left(\frac{H_{ESS}}{L_{min}}\right)^{b}$$

At the model's emergent height, compare $D_{implied}$ against measured basal vessel diameters for the species. If TF24 predicts basal conduits well outside the observed envelope, something upstream is wrong. Same number, different epistemic slot: a falsification test carrying no fitted parameter and introducing no discontinuity.

### 4.3 Elasticities

| Parameter | $\partial \ln R_L / \partial \ln x$ |
|---|---|
| `ks_tip` | −1 (exact) |
| `theta_tip` | +1 (exact) |
| H | $1-\beta$ |
| `L_min` | $+\beta$ |
| b | $\displaystyle \left(2 + \frac{\partial a_\theta}{\partial b}\right)\!\left[\frac{1}{1-\beta} - \ln\frac{H}{L_{min}}\right]$ |

The b-elasticity contains a near-cancellation between the $1/(1-\beta)$ prefactor and $\ln(H/L_{min})$, so it must be evaluated numerically rather than reasoned about — evaluate it before designing the b sweep.

### 4.4 Back-compatibility with the current setup

#### Exact reproduction — the regression gate

Set **b = 0, a_θ = 0, L_min → 0**. Then β = 0 and the closed form collapses to

$$R_L = \frac{\theta_{tip}}{k_{s,tip}}\left(H - L_{min}\right) \;\longrightarrow\; \frac{\theta}{k_s}H$$

bit-identical to current TF24. One coding requirement: implement the integral as

$$R_L = \frac{\theta_{tip}}{k_{s,tip}}\cdot\frac{L_{min}^{\beta}\left(H^{1-\beta}-L_{min}^{1-\beta}\right)}{1-\beta}$$

rather than as $(H/L_{min})^{1-\beta}$, so that `L_min` = 0 is admissible without a division. Guard β = 0 and β = 1 explicitly.

**Correction — do not implement that form.** It fails on two counts. As β → 1 both powers approach 1 and differencing them destroys every significant digit (at 1−β = 1e-12 it retains about four). And its `L_min` = 0 admissibility rests on `pow(H, 1.0) == H`, which is a libm courtesy rather than an IEEE-754 or C99 Annex F guarantee — not something to hang the I7 regression gate on. Implement instead

$$L_{eff} = L_{tip}\cdot\frac{\mathrm{expm1}\!\left((1-\beta)\log(L_{top}/L_{tip})\right)}{1-\beta}$$

which holds full relative precision down to $|1-\beta| \sim 10^{-300}$, needs no separate branch for β > 1 (both numerator and denominator go negative and the quotient stays positive, tending to $L_{tip}/(\beta-1)$ — the saturating regime of §4.1), and takes the logarithmic limit at β = 1 exactly. `L_tip` = 0 is then admissible *only* at β = 0, which is correct: at β > 0 a zero tip length sends $k_s(L)$ to infinity everywhere and is a degenerate configuration to reject loudly, not an edge case to tolerate. The β = 0 branch returns its operand with no arithmetic performed at all, which is what makes I7 bitwise rather than approximately true.

**The path runs to $H\eta_c$, not to $H$.** TF24's existing conductance is `K_s * theta / (height * eta_c)`, where $\eta_c$ = 0.8862 at the default `eta` = 12 is the leaf-area-weighted mean height fraction of the Yokozawa crown. So the representative flow path ends at the mean leaf height. $\eta_c$ scales the **upper limit** of the integral, not the resistance: the two readings differ by the constant $\eta_c^{-\beta}$ (1.05× at β = 0.4), which is H-independent and therefore absorbed wholesale by the reparameterisation below. Every ratio quoted in this section changes slightly as a result — the 9.5× becomes **8.56×** at $H_a$ = 16.6 m.

*Deferred, with the number attached:* the exact leaf-area-weighted $\int f(u)/L_{eff}(uH)\,du$ differs from this mean-path approximation by −0.63% at H = 0.3 m and −0.62% at H = 40 m (β = 0.2, `eta` = 12), shrinking as β grows. Sub-1%, near-constant in H, hence absorbed by the anchor — and the *current* model already carries it at β = 0, so nothing new is introduced. It scales with `eta`, not β, so the check for the mulga growth form (open item 2) is "what `eta` does mulga get".

#### Only one degree of freedom exists

$R_L$ depends on `theta_tip` and `ks_tip` **only through their ratio**. There is exactly one free scalar, so current behaviour can be matched at exactly *one* reference height — never two. Any β > 0 is back-compatible at a point, not across the size range. That is by construction, since changing the height dependence is the object of the exercise.

Matching at $H_a$ requires

$$\frac{\theta_{tip}}{k_{s,tip}} = \frac{\theta_{old}}{k_{s,old}}\cdot\frac{H_a}{G(H_a)}, \qquad \frac{H_a}{G(H_a)} \approx (1-\beta)\left(\frac{H_a}{L_{min}}\right)^{\beta}$$

#### The consistency check this yields

For β = 0.4 (b = 0.2, constant θ), $H_a$ = 20 m, `L_min` = 0.02 m:

| Quantity | Value |
|---|---|
| Required reduction in `ks_tip` vs current bulk $k_s$ | $0.6 \times 1000^{0.4} \approx \mathbf{9.5\times}$ |
| Independent anatomical expectation, $(H_a/L_{min})^{2b}$ | $\approx \mathbf{15.9\times}$ |

Same order, and the gap is the $(1-\beta)$ factor — as expected, since the old uniform $k_s$ behaves like a path-integral average rather than a basal value. **The back-compatible `ks_tip` lands where anatomy says a tip conductivity should sit.**

This is a cheap and genuine check: if a measured twig $k_s$ comes in an order of magnitude away from ~9.5× below the current calibrated value, the old calibration was absorbing something other than path length, and that must be understood before Phase 3.

For β = 0.8 the factor is ~50×, but that case also moves θ, so it is not directly comparable to the anatomical ratio.

#### The second leak, missed by this document

`K_s` is not hydraulics-only either. `make_TF24_hyperpar` derives the entire vulnerability curve from it — $p_{50} = 10^{B_{Hv1} + B_{Hv2}\log_{10}K_s}$ with $B_{Hv2} = -0.2$, and then `c`, `b` and `psi_crit` from $p_{50}$. Dividing `K_s` by 2.98 to get a tip conductivity therefore moves $p_{50}$ from **2.889 to 3.593 MPa**, re-shaping the safety margin of a model whose whole subject is hydraulic limitation. The reparameterisation would have bought a height exponent and silently sold the vulnerability curve.

(Note also that 2.889 MPa is the hyperparameterised $p_{50}$ at `K_s` = 1. `TF24_Pars` defaults `p_50` to 1.85, but `make_TF24_hyperpar` overrides it, so anything quoting 1.85 is quoting the un-hyperparameterised path.)

**Fix: re-anchor the intercept**, $B_{Hv1}: 0.4607063 \rightarrow 0.36591565341924093$, holding $p_{50}$ at 2.8887256653 exactly at the new default. Two consequences to document rather than fix: the roxygen "p50 at `K_s` = 1" now describes a different physical quantity, and `inst/scenarios/scenario_mapping.csv`'s `K_s` High/Low rows must be rescaled to the new baseline or those scenarios silently change meaning.

Whether the Ks–p50 relation should be keyed on tip conductivity at all is a real question — twig segments are what people actually measure — and it needs an explicit closure date. If it was fitted on terminal branch segments, keying it on tip conductivity is *more* correct than before and the re-anchoring should be undone rather than kept.

#### The staged structural/hydraulic θ split — RETRACTED

An earlier version of this section proposed introducing θ(L) for **hydraulics only**, leaving structural θ constant, and unifying them at a later phase. **That staging is retracted and must not be implemented.**

θ has one meaning. It is leaf area per unit sapwood area, and TF24 reads it in both the hydraulic term and the carbon budget — `area_sapwood`, `area_bark`, their growth rates, `mass_sapwood` (hence construction cost, respiration, turnover and NSC capacity), and the hard-coded `dmass_sapwood_darea_leaf` derivative. Splitting it would produce a plant that *conducts* as though θ varied along the stem while being *built and respired* as though it did not: two different plants sharing one trait name. That is not a staging decision, it is an incoherent model, and the fact that it would have been convenient is not a reason.

The consequence is that **θ(L) is a single change, not two**: when it lands, it lands in the hydraulic path *and* in the allometry at the same time. That is harder than the staged version, and it is the honest cost of the compensation mechanism.

Until then `theta_c` is **declared and refused**: `prepare_strategy()` throws on any non-zero value rather than applying a half-model. Widening (`D_c`) carries the height dependence in the meantime, and is unaffected — it enters through $k_s(L)$, which has no structural counterpart to keep in step.

⚠️ When θ(L) is implemented, `dmass_sapwood_darea_leaf` (`src/tf24_strategy.cpp`) is the trap: it hard-codes the closed-form derivative $\rho\eta_c a_{l1}\theta(a_{l2}+1)A^{a_{l2}}$ rather than differentiating `mass_sapwood`, so it must be re-derived by hand and **nothing in the test suite will fail if it is not.**

#### Recommended starting configuration

```
b        = 0.2      # measured, not fitted
a_theta  = 0        # Phase 2a: no theta profile, so `theta` keeps its
                    #   current whole-plant meaning and value unchanged
theta    = unchanged from current  # narrowing to tip definition lands at 2b
L_min    = 0.02     # bound to the ks_tip / theta terminal-segment definition
H_anchor = mature stature (document the value)
ks_tip   = ks_old / 9.5   # then check against measured twig ks
```

Note that `theta` is untouched in Phase 2a by design: with $a_\theta = 0$ there is no profile, so the trait's meaning and numeric value are identical to the current model and no parameter-file migration is needed yet. The ~16× redefinition and the schema bump arrive together at Phase 2b.

Then move `a_theta` from 0 → 2b as a **sweep, not a switch**, so that P2 and P7 are tested rather than assumed.

---

## 5. Gravity — keep as a separate term

$$\Psi_{leaf} = \Psi_{soil} - E\,R_L(H) - \rho g H$$

**Do not lump gravity into an effective resistance.** Three reasons:

1. It multiplies nothing — it exists at E = 0, so lumping breaks the cost function under varying E.
2. The predawn Ψ gradient is pure gravity, and is your cleanest identifiability handle for separating the two terms.
3. Once widening flattens the frictional term to $H^{0.2}$, **gravity carries essentially the whole height signal.** This is not a patch; it is the correct physics, and it makes the height penalty environment-dependent rather than uniform.

Magnitudes: 0.05–0.10 MPa across mulga stature; ≈0.69 MPa at 70 m. One mechanism spanning both endpoints, no tuning.

---

## 6. Integration with the hydraulic cost component

Scope: the cost term only, per the variant ProfitMax formulation in use.

$$\mathcal{C}(E) = 1 - \frac{k(\Psi_{leaf})}{k_{max}}$$

### 6.1 The critical structural point

**The cost term is normalised, so it is invariant to uniform rescaling of conductance.** $k/k_{max}$ is dimensionless and depends only on where $\Psi_{leaf}$ sits on the *shape* of the vulnerability curve. Consequences:

- Changing `ks_tip`, or changing the height exponent from 1 to 0.2, **rescales flux but barely moves the cost term.** Any residual effect enters only through curvature of $A(g_c)$ on the gain side.
- Therefore **all first-order height sensitivity in the cost term must come from ρgH**, which translates the operating point along the Ψ axis rather than rescaling it.

$$\Delta E_{crit} \approx \int_{\Psi_s - \rho g H}^{\Psi_s} k\,d\Psi \approx k_{max}\,\rho g H$$

Gravity truncates the **high-conductance plateau first** — the most valuable part of the pressure envelope, at the maximum possible rate.

### 6.2 Lumped vs segmented vulnerability

Recommend **lumped, single vulnerability curve, anchored to tip conduits.** Rationale: multiple $P_{50}$ values are not identifiable from flux data; and Ψ is most negative at the tips, so that is where cost is realised.

Known cost of lumping, to log rather than fix: axial $P_{50}$ gradients are real (in Douglas-fir, roughly −4.7 MPa near the apex to −3.3 MPa at the base), so a lumped curve will misplace where failure initiates. Defer.

### 6.3 Optional couplings to sweep, not fit

- $P_{50}(D)$: conduit-level safety–efficiency coupling is empirically weak, so a null coupling is defensible as the default. Sweep to check whether conclusions depend on it.

### 6.4 Calibratability trap — do not do this

Setting `ks_tip` (or $k_{max}$) at each height by requiring an optimal $c_i:c_a$ ratio absorbs the entire hydraulic limitation into the parameterisation. The model then fits well and predicts nothing about height, because the mechanism under test has been used as the calibration target. If the coordination hypothesis is wanted, calibrate at **one** reference height and let $H^{1-\beta}$ propagate.

---

## 7. Predicted emergent behaviours

Each of these is a claim to demonstrate over sweeps, not to fit.

**P1 — The marginal cost of height rises with aridity.** ρgH is a fixed absolute offset but a variable *fraction* of the available pressure envelope: small when Ψ_crit is far away, large when drought has compressed the envelope. This gives Paper G an aridity-dependent height penalty that is derived rather than tuned.

**P2 — Emergent ontogenetic Huber trend.** $1/\theta$ should increase with H as $H^{a_\theta}$. With $a_\theta = 2b$ this is a prediction testable against McDowell et al.'s Douglas-fir series, not an input. Check against the derived output `theta_base(H)` (§1), which is the model quantity directly comparable to whole-tree A_l:A_s measurements.

**P3 — Emergent *negative* Hv–H_max relation across strategies.** Within species Huber value rises with height; across species, taller-statured species have *lower* Huber values, attributed to sapwood carbon cost versus leaf gain rather than to hydraulics. These pull in opposite directions. TF24 already contains the relevant carbon economics — if it *derives* the negative across-strategy relation while the ontogenetic profile runs the other way, that is a Paper A′ result. If it is imposed, the result has been spent to buy a parameter.

**P4 — struck.** An earlier draft predicted candidate branching driven by a widening cap. This was circular: the cap was introduced in §2.1 and the prediction then built on it, and the derivative discontinuity a cap creates generates singular strategies at the kink as an artefact of non-smoothness rather than as biology. Removed rather than flagged. If branching in height strategy does emerge from the smooth formulation, it needs a named mechanism.

**P5 — ESS height may become unbounded.** Linear-R was probably doing hidden work bounding optimal height. Weakening it ~4-14× may reveal that nothing else bounds it. If so, that is diagnostic, not a failure: it tells you linear-R was standing in for an unspecified cost (sapwood respiration scaling with sapwood volume, mechanical investment, or self-shading) that now needs naming. **Better found deliberately in a sweep than discovered as a runaway.**

**P6 — Establishment shifts.** Reparameterisation changes seedling resistance materially depending on anchor height. Recruitment and hence emergent size distributions will move. Must be checked before any demographic conclusion.

**P7 — Invariance at the arid endpoint (the key test).** At 5–10 m, ρgH ≈ 0.05–0.10 MPa and the $H^{0.2}$ term is nearly flat. **The mechanism ranking at Alice Mulga / Ti Tree East should therefore be invariant to `b` and `a_theta` across their full plausible ranges.** If it is not, something is wrong upstream. This is a stronger protection argument than reporting an improved fit.

**P8 — Finite height limit at the mesic endpoint without imposition.** With the cap removed, gravity truncating the pressure envelope must do this work, together with sapwood respiration scaling with sapwood volume and the carbon cost of widening itself. If it does not, see P5 — the answer is to name the missing cost, **not** to reinstate a diameter cap.

---

## 8. Empirical support, by link

| Link | Evidence | Status |
|---|---|---|
| b ≈ 0.2 | Conserved across terrestrial vascular plants; optimality models jointly minimising resistance, construction cost and embolism risk predict roughly that value. Widened pipe model (Koçillari et al. 2021). Within-stem range ~0.1–0.3. | **Strong** |
| Widening in very tall trees | *Sequoia*, *Sequoiadendron*, *E. regnans* at 86–105 m: power-function decline in widening per unit path length, consistent with hydraulic compensation but bounded by maximum conduit diameter. | Observation strong; **mechanism unknown, so not used as a model constraint** (§2.1, §4.2) |
| Over-compensation in *E. regnans* | Measured basal $D_h$ 220–250 µm vs 160–175 µm predicted by metabolic scaling theory — implies optimisation for hydraulic efficiency over carbon cost. | Moderate; single study |
| $n_A \propto D^{-2}$ (packing limit) | Strong across species (Sperry, Hacke & Pittermann 2008). Less directly tested *along an axis within a tree*. | Strong interspecific, **assumed within-axis** |
| Lumen fraction conserved along axis | Interspecific conservation documented; basal wood density gradients could violate it, which would worsen compensation. | **Assumed, untested** |
| Path length (not cambial age) drives conduit size | Bicego et al. (2025) coppice manipulation: median conduit area at breast height fell in all 7 sprouting trees (ratios 0.93–0.56, all P < 0.01) while 2 unsprouted controls were unchanged (1.00, 1.04). Runs the experiment *backwards*, dissociating path length from age. Fajardo et al. (2020): stem length not climate controls vessel diameter. | **Strong causal, n = 9, one site** |
| Huber value rises with H *within* species | McDowell et al. (2002): 15 Douglas-fir, 13–62 m, aged 20–450 y; A_l:A_s declined substantially (P = 0.02), extended across nine species. Confirmed as the within-species pattern by Mencuccini et al. (2019). | Strong |
| Huber value falls with H_max *across* species | Mencuccini et al. (2019): negative isometric relation, attributed to sapwood carbon costs vs leaf gains. | Strong, and opposite in sign to the above |
| $k_L$ declines with height | Douglas-fir: 44% decrease from 15 → 32 m, then only 6% further to 60 m — strongly saturating. Zaehle (2005): compensation incomplete; *Betula occidentalis* >2-fold decline. Ryan et al. (2006), 51 studies: reductions common but not universal, complete compensation rare. Counter-case: *E. saligna* at 7 vs 26 m fully compensated. | Strong but heterogeneous |
| Saturating rather than power-law form | Zhao et al. (2026), 141 poplar saplings: leaf-area-basis resistance increased with size; asymptotic models outperformed power laws. | Moderate, recent |
| **Direct hydraulic measurement of R vs path length** | Mostly anatomy + Poiseuille integration, not measured flow. Segment-level measurements exist (Petit et al. 2008 sycamore; Pfautsch et al. 2018 *E. grandis*; AoB PLANTS 2025 tip-to-base). **No within-individual longitudinal whole-plant K measurements exist.** | **Weak — the main gap** |
| Ultra-widening permeability | Anfodillo & Olson (2024): stretching argument implies a steeper requirement than widening alone delivers; hypothesised pit membrane thinning. | **Unresolved — log as residual, do not absorb into `ks_tip`** |
| Gravity limits height | Koch et al. (2004), redwood. | Strong |
| Cost formulation without fitted parameters | Sperry et al. (2017). Sicangco et al. (2026): alternative implementations tested against *E. parramattensis* whole-tree chamber heatwave data. | Strong; read Sicangco before fixing the cost form |

**Epistemic note for the ledger:** the data here are conduit diameters. Height-invariant conductance is a *phenomenon inferred through a model*, not an observed quantity. The Bogen–Woodward distinction should be explicit wherever the $H^{0.2}$ result is cited.

---

## 9. Why capacitance can be deferred

Given the target timescale is environmental change over weeks to years, capacitance can be omitted, and this should be stated as an explicit scope assumption rather than a silent omission.

**Argument:**

1. **It is a closed diurnal cycle.** Storage discharges by day and recharges overnight, so net change over 24 h is zero. At daily resolution there is nothing left to represent. Meinzer's own scope condition is the operative one: storage produces lags between transpiration and basal sap flux but is not expected to affect the size dependence of total daily transport unless overnight recharge is incomplete.

2. **It cannot buffer gravity at all.** ρgH is a standing offset present at zero flow — measurable predawn. A reservoir that is full and static cannot offset a term that exists when it is full and static. Capacitance shaves the diurnal peak of the *frictional* term only, and the frictional term is the one widening already handles.

3. **The buffer is self-limiting.** Stomatal closure from the early-morning maximum limits tension and thereby prevents further reservoir discharge, which is why some studies find storage contributing a roughly constant ~10% fraction regardless of size.

4. **A larger error is already present.** Collapsing the diurnal cycle to daily means introduces a Jensen error across the concave $A(g_c)$ curve that almost certainly exceeds the capacitance term. Optimising capacitance while using daily-mean VPD would be fixing the wrong term. **Priority: get the diurnal closure right first** ([#618](https://github.com/traitecoevo/plant/issues/618)).

5. **Under well-watered conditions the effect is small anyway** — a 30-fold sweep in stem capacitance moved GPP by only ~0.2 g C m⁻² day⁻¹ at high soil moisture.

**The one exception, and it is not diurnal.** Incomplete overnight recharge across a dry-down is a genuinely multi-day mechanism: a slow reservoir depleting cumulatively, setting the lag between soil water decline and stomatal closure. Timescale days to a fortnight — exactly the Alice Mulga pulse timescale. If any storage term earns its place in Paper G it is this one, and it is a *different* parameterisation (a slow state with a recharge rate, not a capacitor with τ = RC). Under low soil moisture the same capacitance sweep moved GPP by ~0.67 g C m⁻² day⁻¹, roughly three times the wet-condition sensitivity — a state-dependent sensitivity that no single calibrated conductance can reproduce. Treat as a scoped sweep, not a default component.

**Scope assumption to log:** complete overnight recharge — assumed, untested at Alice Mulga and Ti Tree East.

### 9.1 Required companion work — moved out

The diurnal closure — that 𝒞 needs the daily *minimum* Ψ_leaf while daily A is an *integral* over a concave $A(g_c)$, so a single daily Ψ serves neither — is **not part of this plan**. It is temporal aggregation, not stem architecture, and it does not depend on the work here.

Tracked separately as [#618](https://github.com/traitecoevo/plant/issues/618). Note that argument 4 above leans on it: the claim that capacitance can be deferred rests on the Jensen error being the larger term, which #618 is what would establish.

---

## 10. Implementation sequence

| Phase | Work | Gate |
|---|---|---|
| 0 | Naming: `b`→`D_c`, `a_theta`→`theta_c`, `L_min`→`L_tip`; `K_s` keeps its name. Cost-θ rename is a no-op (§1) | **done** |
| 1 | Implement $k_s(L)$, $n_A(L)$, $\theta(L)$ profiles; closed-form $L_{eff}$ with `theta_c` as switch, all three parameters defaulting to 0 | reproduces linear case **bit-identically**; `scientific_version` NOT bumped |
| 2a | Turn widening on (`D_c` = 0.2, `L_tip` = 0.02) and reparameterise `K_s` at `H_anchor` = 1 m; `theta_c` stays 0 and is refused; re-anchor `B_Hv1` to hold $p_{50}$ | resistance unchanged at $H_a$; carbon budget untouched because θ has no profile yet; vulnerability curve unchanged; `scientific_version` → 10 with measured deltas |
| 2b | Implement θ(L) **everywhere at once** — hydraulic path *and* allometry (`area_sapwood`, `area_bark`, `mass_sapwood`, `dmass_sapwood_darea_leaf`, and their rate forms); narrow `theta` to the tip definition; **build** the parameter-file schema gate; expose `theta_base(H)` as derived output | `theta_c` guard lifted; carbon budget and hydraulics use the same θ; old-format parameter files rejected (I11) |
| 3 | Add ρgH as an explicit separate term | predawn Ψ gradient reproduced |
| 4 | Diagnostic ceiling check (§4.2) | $D_{implied}$ at emergent height compared against measured basal diameters |
| 5 | Sweeps and invariance tests | §11 criteria met |

### 10.1 Sweeps

| Axis | Range | Purpose |
|---|---|---|
| `b` | 0.10 – 0.30 | observed within-stem range |
| `a_theta` | 0 – 0.6 | constant θ → saturating |
| `H_anchor` | seedling → mature | establishment sensitivity (P6) |
| ρg | on / off | isolate gravity's contribution (P1) |
| $P_{50}(D)$ coupling | null / weak | safety–efficiency dependence |

---

## 11. Invariance criteria

- **I1** Mechanism ranking at Alice Mulga / Ti Tree East is invariant to (`b`, `a_theta`) over their full plausible ranges.
- **I2** Emergent ontogenetic Huber trend falls within the confidence interval of McDowell et al.'s data *without* having been fitted to it.
- **I3** ESS H_max is finite with gravity on and $a_\theta = 2b$; if it is not, the binding cost is identified and named rather than restored by re-steepening R.
- **I4** No parameter is calibrated to a target downstream of the mechanism under test (specifically: `ks_tip` is not set via a $c_i:c_a$ criterion at more than one height).
- **I5** `ks_tip`, `theta` and `L_min` all originate from a single terminal-segment definition.
- **I6** Any conductance shortfall relative to observation is recorded as an explicit ultra-widening residual, not absorbed into `ks_tip`.
- **I7** Regression: at b = 0, a_θ = 0, `L_min` = 0 the implementation reproduces current TF24 output exactly (§4.4).
- **I8** The back-compatible `ks_tip` at $H_a$ agrees to within an order of magnitude with a measured twig $k_s$, and the discrepancy is accounted for rather than tuned away.
- **I9** $D_{implied}$ at the emergent height falls within the observed basal-diameter envelope for the focal species (§4.2). Violation is a falsification signal, not a licence to add a cap.
- **I10** $R_L(H)$ is continuously differentiable over the whole height range. Any kink introduced later must be traced to a named mechanism before it is allowed to influence ESS height.
- **I11** No parameter file written under the old `theta` definition can be loaded once $a_\theta > 0$: the schema version gate rejects rather than reinterprets it. ⚠️ The gate does not exist; Phase 2b must build it (§1).
- **I12** $p_{50}$ at the default parameters is unchanged by the `K_s` reparameterisation. Whatever the reparameterisation buys, it must not also move the vulnerability curve (§4.4, the second leak).
- **I13** θ has exactly one meaning. The hydraulic term and the carbon budget read the same θ(L); no configuration may profile one without the other (§4.4).

---

## 12. Open items

1. `b` elasticity contains a near-cancellation — evaluate numerically before designing the sweep.
2. Crown-thickness assumption in the §4 closed form: valid when the crown is thin relative to H. Check for the mulga growth form, where it may not be.
3. Root and rhizosphere conductance is untouched here and is plausibly what actually limits the arid endpoint. The height-hydraulics apparatus above is near-irrelevant at 5–10 m; do not let it absorb work that belongs to the soil–root term.
4. Obtain measured basal vessel diameters for the focal species to make the §4.2 diagnostic ceiling operational. Without them the check cannot be run.
5. **Closure date for the Phase 2a structural/hydraulic θ inconsistency.** Set one now (§4.4). This is the item most likely to become permanent by neglect.
6. `theta` and `ks_tip` are identifiable only as a ratio, so back-compatibility holds at one height only. Decide and document `H_anchor` before Phase 2a, since it is not recoverable afterwards from model output alone.
7. Obtain a measured twig $k_s$ for at least one mulga and one tall-eucalypt species to close I8. Without it, §4.4's consistency check is unusable.
