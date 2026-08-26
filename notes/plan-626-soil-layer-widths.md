# Variable soil layer widths (#626)

Status: **implemented**, in plant (`feature/soil-layer-widths-626`) and phylloptim 0.9.0 (`feature/per-layer-dz-626`). Base branch is `feature/events-522` (PR #632), not `develop`.

Soil evaporation needs a very thin layer at the top of the profile, which a layer *count* cannot express. [#626](https://github.com/traitecoevo/plant/issues/626) asks for variable widths and prescribes the API: **input the width of each layer, not the depth**, because a "depth" is ambiguous between a top, a midpoint and a bottom.

## The framing question: layers or nodes?

It was put that plant "doesn't really have layers, but nodes at depth". Half true, and the half that is false is the structural half.

`TF24_Environment` was already a finite-volume **layer** discretisation carrying all three geometric descriptions, and each is read by exactly the part of the model that needs it:

| member | meaning | read by |
|---|---|---|
| `dz[i]` | layer **width** — already a vector | the water balance; the pulse capacity cap |
| `z[i]` | cumulative depth to the **bottom** of layer `i`; `z.back() == depth` | root-mass distribution, `Q(z)` differenced across a layer |
| `z_mid[i]` | layer **midpoint** | the leaf's gravitational head |

So the values were not "none of top/midpoint/bottom" — `z[i]` is unambiguously the bottom, and both packages treated it that way consistently. The state itself is θ, intensive, so it needs no width at all.

The confusion was **nominal**: the setter was `set_soil_number_of_depths`, the R-visible vector was `soil_depth`, `delta_z` was commented "distance between layers" (a nodal concept) though used as a width, and R could not ask for the geometry at all. **The fix is definitional** — width is the only geometry a caller states, boundaries and midpoints are strictly derived, and nothing is left to be a top, a midpoint or a bottom.

## Decisions, and what constrains future edits

### 1. The inter-layer cascade needs no change. Do not "generalise" it.

`water_flux[i] = K(theta_i)`, with layer `i+1` receiving it, is **free drainage**: Darcy's `q = -K(dpsi/dz + 1)` with the matric gradient dropped. Because that term is absent, **no internode distance appears**, so there is nothing to rescale for unequal layers. It is evaluated at the donor cell, which is correct upwinding. Mass conservation `sum(dz[i]*dtheta_i/dt) = infiltration - drainage - uptake` is exact at any width distribution, and is now pinned on four graded profiles in `test-tf24-water-budget.R`.

⚠️ **If anyone adds a capillary or diffusive flux — and soil evaporation will make them want to — that term's distance is `z_mid[i+1] - z_mid[i]`, NOT `dz[i]`**, and the interfacial conductivity must be a distance-weighted harmonic or geometric mean of the two cells. The two distances coincide only for equal layers, so a uniform-profile budget test passes on the mistake and the graded one does not. Recorded at the `compute_rates` comment block as well as here.

### 2. Two floating-point constraints, both load-bearing

Measured over 180 `(depth, n)` pairs:

| route to a uniform width vector | disagrees with today's values |
|---|---|
| pass the environment's own `dz` | **0 / 180** |
| difference the profile, `diff(c(0, z))` | 122 / 180 |
| rebuild `z` as `cumsum(dz)` | 103 / 180 |

Therefore:

- **`TF24_Strategy` passes `environment.dz` and must never derive widths from `z`.** `phylloptim::layer_thicknesses()` is for a caller holding only a profile; from plant it is actively wrong. Stated at the call site and in `agents.md`.
- **`set_soil_number_of_depths` keeps `z[i] = (i+1)*delta_z` and must not be unified with `set_soil_layer_widths`' running sum.** `z` feeds `Q()` and the gravitational head, so unifying them re-baselines every TF24 number for no gain. The `test-environment-TF24.R` geometry test compares the two paths with `expect_equal`, not `expect_identical`, and says why.

### 3. There was a real 3.7× bug behind this, now fixed upstream

`phylloptim::layer_thickness()` returned the scalar `z.back()/n` and `root_network_from_carbon` squared it. Correct only for equal layers; what fails otherwise is **discretisation invariance** — slicing a column is a numerical choice, so nothing about the plant may depend on it. With root density uniform over depth, layer-integrated carbon goes as `dz[i]`, giving a total `sum(r_R_V) = 3*beta_R_V*D^2/C` for any slicing. Measured, 1.5 m carrying 20 kg C m⁻² leaf:

| profile | scalar `dz` | per-layer `dz[i]` | ratio |
|---|---|---|---|
| 5 or 3 equal layers | 3172.5 | 3172.5 | 1.000 |
| 5 cm surface layer | 6059.5 | 3172.5 | 1.910 |
| 2 cm surface layer | 11688.4 | 3172.5 | **3.684** |

At plant's own operating point (θ = 0.25, height 5 m), measured by building against each:

| profile | `opt_root_psi` (MPa) scalar → per-layer | `assimilation` scalar → per-layer |
|---|---|---|
| uniform, 5 layers | 1.661740 → 1.661740 | 15.464934 → 15.464934 |
| 2 cm surface layer | 1.747280 → 0.598266 | 15.338550 → 17.251842 |
| thick over thin | 1.654759 → 4.211291 | 15.478450 → 11.179142 |

A scalar thickness lands every graded profile within 6% of the uniform value and **on the wrong side of it**, so `test-tf24-root-pars.R` diffs the *sign*, not a tolerance.

## What is deliberately NOT in scope

Soil evaporation. #626 is the enabling geometry change only. The follow-up owns four things, none of them geometry:

1. **Re-key the infiltration switch off a fixed reference depth.** It reads layer 0 alone through `(theta_0/theta_sat)^b_infil` with `b_infil = 8`, and `a_infil`/`b_infil` were calibrated against a 30 cm surface layer. Measured in the R replica (`notes/scripts/soil_cascade/`), 20 events × 50 mm: annual runoff fraction **0.134 → 0.204** from thinning layer 0 alone, no parameter changed. Inert under constant rain — the steady state contains no `dz` — so the whole effect is in the event transient. **This is the item most likely to be missed.**
2. **Step control for a thin wet layer** — see below; the concrete case for the IMEX/substepping option deferred in `plan-tf24-soil-redistribution.md`, interacting with #599 and #529.
3. **Upward capillary flux**, which is the term that introduces the internode distance (§1).
4. **The evaporation sink and its driver.** `wind_speed` and the Penman-Monteith machinery from #523 already exist in the environment.

## The stiffness limit, which this change does not fix

Shipped **permissive**: widths are validated for finiteness and positivity and nothing else. The numbers, so they are on the record:

| min `dz` | diagonal `K'(theta)/dz` at 0.99·θ_sat | stable step `2/lambda` |
|---|---|---|
| 0.30 m (default) | 1.8e4 yr⁻¹ | 1.1e-4 yr |
| 0.10 m | 5.3e4 | 3.8e-5 |
| 0.02 m | 2.6e5 | 7.6e-6 yr ≈ 4 min |
| 0.01 m | 5.3e5 | 3.8e-6 |

The analytic scaling reproduces the 1.77e4 / 5.30e4 measured in `plan-tf24-soil-redistribution.md`. Two consequences:

- It only bites when the layer is **wet**: at dz = 0.02 the diagonal is 2.6e5 at 0.99·θ_sat but 8.5 at 0.5·θ_sat. So the hazard is the post-rain transient — exactly when an evaporation layer matters.
- `Control()$ode_step_size_min` is 1e-6 yr, which the stable step reaches at **dz ≈ 2.6 mm**. Below that the controller cannot satisfy stability at all. That is a hard floor on layer thickness in this model.

`soil_widths_graded()`'s roxygen carries both this and the runoff caveat, since that is where a user deliberately asks for a thin layer.

## Verification performed

- phylloptim: full suite green; all three golden baselines (576 operating points, 5184 psi_stem optima, 544 primitives) **bit-identical**; `gradient_golden.tsv` bit-exact; both bench programs build; the CI consumer program (a heredoc inside `.github/workflows/cpp-tests.yml`, which neither `make` nor `cmake` covers locally) extracted, compiled and run.
- plant: `make clean && make RcppR6 && make compile`; full `test_dir` green; `test-model-version.R` (the scientific-surface hash) green with `NOT_CRAN=true`.
- **Bit-exactness on the default geometry, at SCM level:** `run_scm` on the default uniform environment gives `11.143684537309916` under both the scalar and the per-layer build. The same graded profile gives 9.2796204566000107 (scalar) against 12.981712222339016 (per-layer).
- ⚠️ **The scenario gateway is red, and it was red before this change.** Every scenario's `offspring_production` differs from `tests/testthat/test_data/scenario_baseline.rds` by 19–41× on S03–S06 — the range #590's density-coordinate change is documented to produce — so the recorded baseline predates it. Do **not** `make bless-scenarios` as part of this change: the bit-identical SCM result above is what says #626 did not move it. Re-blessing is its own decision, on its own branch.

## Release order

phylloptim lands first. Its pin is `==`, and plant compiles against the *installed* headers, so:

1. merge phylloptim 0.9.0, record the merge sha;
2. install it locally **before** touching plant's pin (bumping the pin over a stale install fails in a way that reads like a broken build — it did here);
3. plant `DESCRIPTION`: `LinkingTo: ... phylloptim (== 0.9.0)` and `Remotes: traitecoevo/phylloptim@<merge-sha>`;
4. `make clean && make rebuild` — clean is mandatory.

`scientific_version` is **not** bumped: no default output moves.

⚠️ **Stacking.** This is based on `feature/events-522` (PR #632), which will squash-merge. After it lands, `git rebase --onto develop feature/events-522 feature/soil-layer-widths-626` and check the changed-file count is ~15, not ~55, before merging — otherwise the child re-lands #632's whole 2800-line diff.
