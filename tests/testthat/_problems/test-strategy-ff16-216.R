# Extracted from test-strategy-ff16.R:216

# test -------------------------------------------------------------------------
p0 <- scm_base_parameters("FF16")
env <- Environment("FF16")
ctrl <- Control()
p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
out <- run_scm(p1, env, ctrl)
expect_equal(out$offspring_production, 16.88946, tolerance=1e-4)
expect_equal(out$ode_times[c(10, 100)], c(0.000070, 4.216055), tolerance=1e-5)
