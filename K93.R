devtools::load_all()

pl <- K93_Plant()
pl$state("height")

pl$set_state("height", 10)
env <- K93_fixed_environment(1.0)
pl$compute_rates(env)

devtools::load_all()
p0 <- scm_base_parameters("K93")
p1 <- expand_parameters(trait_matrix(0.059, "b_0"), p0, K93_hyperpar, FALSE)
p1$seed_rain <- 20
out <- run_scm(p1)


# one species

p0 <- scm_base_parameters("FF16")
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16_hyperpar, FALSE)

p1$seed_rain <- 20
out <- run_scm(p1)
expect_equal(out$seed_rains, 16.88946, tolerance=1e-5)
expect_equal( out$ode_times[c(10, 100)], c(0.000070, 4.216055), tolerance=1e-5)
