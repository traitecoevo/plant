## Rscript probe_inv.R <label>: test-scm.R's invader-vs-stand check, with the gap.
suppressPackageStartupMessages(library(plant))
label <- commandArgs(TRUE)[1]
p0 <- scm_base_parameters("TF24", "TF24_Env")
p0$max_patch_lifetime <- 3
p0$patch_type <- "fixed"
p <- add_strategies(p0, trait_matrix(0.1978791, "lma"))
env <- function() {
  e <- Environment("TF24")
  t <- seq(0, 3, length.out = 601)
  e$extrinsic_drivers_set_variable("rainfall", t, 0.25 * (1 + 0.95 * sin(2 * pi * t)))
  e
}
ctrl <- Control(node_density_in_birth_date = TRUE, ode_tol_abs = 1e-8,
                ode_soil_substep = TRUE)
s <- run_scm(p, env(), ctrl, keep_states = TRUE)
ratios <- s$net_reproduction_ratios
s$run_mutant(p)
inv <- s$net_reproduction_ratios
cat(sprintf("%-14s stand %.17g invader %.17g rel %.3e\n", label, ratios, inv,
            abs(inv - ratios) / abs(ratios)))
