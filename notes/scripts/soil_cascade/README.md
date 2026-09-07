# Soil-cascade experiments (Stage A design pass, #522 / #599)

Standalone R replica of TF24's soil-water cascade — no `plant` build required. Backs the
findings in [`../../plan-tf24-soil-redistribution.md`](../../plan-tf24-soil-redistribution.md).

Run from this directory:

```sh
Rscript exp1.R    # vector field, Jacobian structure, spectrum, candidate schemes
Rscript exp2.R    # odelia's exact controller; the NaN-acceptance bug; pulse capacity
Rscript exp4.R    # the decisive one: inherited step size x throw/NaN leaf handoff
```

| file | contents |
|---|---|
| `soil_cascade.R` | mirrors `tf24_environment.h` `compute_rates()` and `soil_K_from_soil_theta()`, plus the three candidate redistribution schemes (`baseline`, `recv`, `donor`), a numerical Jacobian, and a Cash-Karp integrator |
| `exp2_ctl.R` | odelia's `OdeControl` (`errlevel`, `adjust_step_size`) transcribed exactly, with plant's defaults from `src/control.cpp`, and a stage evaluator that recomputes uptake at every RK stage the way `Patch` does |
| `exp1.R` | Q1–Q5: is θ_sat a barrier; reproducing the overshoot; Jacobian structure and eigenvalues per scheme |
| `exp2.R` | Q6–Q9: does odelia's controller reject the overshoot; the NaN path; pulse-jump capacity |
| `exp4.R` | Q12–Q13: entering a leg with an inherited large step, crossed with a leaf that throws vs returns NaN, with and without the controller fix |

Parameters are TF24 defaults (θ_sat = 0.428, K_sat = 163.0411, n_psi = 6.57, a_infil = 1,
b_infil = 8, depth = 1.5 m / 5 layers) and plant's ODE control defaults (`tol_rel` =
`tol_abs` = 1e-4, `a_y` = 1, `a_dydt` = 0, `h_min` = `h_init` = 1e-6, `h_max` = 5). If any of
those move in the package, update `soil_cascade.R` / `exp2_ctl.R` to match before rerunning.
