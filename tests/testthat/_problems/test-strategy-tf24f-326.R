# Extracted from test-strategy-tf24f.R:326

# test -------------------------------------------------------------------------
mk <- function(strat) {
    p0 <- scm_base_parameters(strat)
    p0$max_patch_lifetime <- 5
    add_strategies(p0, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                   hyperpar = get(paste0(strat, "_hyperpar")), birth_rate = list(20))
  }
set_k_acclim <- function(p, k) {
    s <- p$strategies[[1]]
    s$k_acclim <- k
    p$strategies[[1]] <- s
    p
  }
expect_true("opt_root_psi_state" %in% TF24f_Individual()$ode_names)
pf <- mk("TF24f")
slow <- run_scm(set_k_acclim(pf, 0.1), Environment("TF24f"), Control())$offspring_production
fast <- run_scm(set_k_acclim(pf, 10),  Environment("TF24f"), Control())$offspring_production
expect_length(fast, 1)
expect_true(is.finite(fast) && fast > 0)
expect_gt(abs(fast - slow) / slow, 0.05)
