## Rscript probe_heap.R <label>: the same run repeated in one process, the heap
## disturbed between repeats. A value that moves depends on memory layout.
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
keep <- list()
for (i in 1:6) {
  s <- run_scm(p, env(), ctrl, keep_states = TRUE)
  r0 <- s$net_reproduction_ratios
  st <- s$patch$ode_state
  s$run_mutant(p)
  r1 <- s$net_reproduction_ratios
  cat(sprintf("%s rep %d stand %.17g invader %.17g rel %.3e state_sum %.17g\n",
              label, i, r0, r1, abs(r1 - r0) / abs(r0), sum(st)))
  ## Disturb the heap: hold some odd-sized blocks, and run another model.
  keep[[i]] <- lapply(seq_len(50 * i), function(k) numeric(k * 7 + i))
  f <- add_strategies(scm_base_parameters("FF16"), trait_matrix(0.0825 + i / 1000, "lma"),
                      birth_rate = 1)
  invisible(run_scm(f))
}
