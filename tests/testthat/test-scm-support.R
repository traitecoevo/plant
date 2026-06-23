context("SCM support")


test_that("collect", {
  
  env <- Environment("FF16")
  ctrl <- Control()
  p0 <- scm_base_parameters("FF16")
  p0$disturbance_mean_interval <- 30.0
  p1 <- add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)

  expect_silent(res <- run_scm(p1, env, ctrl))

  expect_contains(
    names(res), c("clone", "collect", "compute_competition_effect_error_by_node_for_species_i", "history", "initialize", "net_reproduction_ratio_errors", "net_reproduction_ratio_for_species", "net_reproduction_ratios", "node_schedule", "ode_times", "offspring_production", "parameters", "patch", "reset", "run", "run_mutant", "set_node_schedule_times", "time")
  )

})

test_that("expand_parameters & mutant_parameters", {
  
  hyperpar <- make_FF16_hyperpar()
  p0 <- scm_base_parameters("FF16")

  expect_equal(p0$strategies |> length(), 0)

  p1 <- add_strategies(p0, trait_matrix(0.1, "lma"), birth_rate = 1.0)
  
  expect_equal(p1$strategies |> length(), 1)
  expect_equal(p1$strategies[[1]]$pars$lma, 0.1)

  p1$max_patch_lifetime <- 100
  expect_silent(p2 <- add_strategies(p1, trait_matrix(0.2, "lma"), birth_rate = 1.0))
  expect_equal(p2$max_patch_lifetime, p1$max_patch_lifetime)

  expect_equal(p2$strategies |> length(), 2)
  expect_equal(p2$strategies[[1]]$pars$lma, 0.1)
  expect_equal(p2$strategies[[2]]$pars$lma, 0.2)

  expect_silent(p3 <- add_strategies(p1, trait_matrix(0.3, "lma"), birth_rate = 1.0, keep_existing = FALSE))

  expect_equal(p3$strategies |> length(), 1)
  expect_equal(p3$strategies[[1]]$pars$lma, 0.3)

  expect_silent(p4 <- add_mutant(p1, trait_matrix(0.3, "lma"), birth_rate = 1.0))

  expect_equal(p3, p4)

})

test_that("collect_auxiliary_variables", {
  
  env <- Environment("FF16")
  ctrl <- Control()
  p0 <- scm_base_parameters("FF16")
  p0$disturbance_mean_interval <- 30.0
  # two species
  p1 <- add_strategies(p0, trait_matrix(0.082, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(11.99177))

  results <- run_scm(p1, env, ctrl, collect = TRUE)
  
  # check columns,should contain auxillary variables
  expect_equal(ncol(results$species), 16)
  expect_contains(names(results$species), c("competition_effect", "height_inverse", "net_mass_production_dt"))
})
