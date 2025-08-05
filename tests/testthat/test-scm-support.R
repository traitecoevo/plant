context("SCM support")


test_that("collect", {
  
  env <- Environment("FF16")
  ctrl <- scm_base_control()
  p0 <- scm_base_parameters("FF16")
  p0$disturbance_mean_interval <- 30.0
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, birth_rate_list = 1.0)

  expect_silent(res <- run_scm(p1, env, ctrl))

  expect_contains(
    names(res), c("clone", "collect", "competition_effect_error", "complete", "history", "initialize", "net_reproduction_ratio_errors", "net_reproduction_ratio_for_species", "net_reproduction_ratios", "node_schedule", "ode_times", "offspring_production", "parameters", "patch", "reset", "run", "run_mutant", "run_next", "set_node_schedule_times", "time", "use_ode_times")
  )

})

test_that("expand_parameters & mutant_parameters", {
  
  hyperpar <- make_FF16_hyperpar()
  p0 <- scm_base_parameters("FF16")

  expect_equal(p0$strategies |> length(), 0)

  p1 <- expand_parameters(trait_matrix(0.1, "lma"), p0, birth_rate_list = 1.0)
  
  expect_equal(p1$strategies |> length(), 1)
  expect_equal(p1$strategies[[1]]$lma, 0.1)

  p1$max_patch_lifetime <- 100
  expect_silent(p2 <- expand_parameters(trait_matrix(0.2, "lma"), p1, birth_rate_list = 1.0))
  expect_equal(p2$max_patch_lifetime, p1$max_patch_lifetime)

  expect_equal(p2$strategies |> length(), 2)
  expect_equal(p2$strategies[[1]]$lma, 0.1)
  expect_equal(p2$strategies[[2]]$lma, 0.2)

  expect_silent(p3 <- expand_parameters(trait_matrix(0.3, "lma"), p1, birth_rate_list = 1.0, keep_existing_strategies = FALSE))

  expect_equal(p3$strategies |> length(), 1)
  expect_equal(p3$strategies[[1]]$lma, 0.3)

  expect_silent(p4 <- mutant_parameters(trait_matrix(0.3, "lma"), p1, birth_rate_list = 1.0))

  expect_equal(p3, p4)

})

test_that("collect_auxiliary_variables", {
  
  env <- Environment("FF16")
  ctrl <- scm_base_control()
  p0 <- scm_base_parameters("FF16")
  p0$disturbance_mean_interval <- 30.0
  # two species
  p2 <- expand_parameters(trait_matrix(c(0.0825, 0.2625), "lma"), p0, FF16_hyperpar,
    birth_rate_list = list(11.99177, 16.51006)
  )

  # Compare results to a reference of prior behavior. 
  # Reference was caulcated from commit 4ca3a9be on develop, before simplifying interface for scm and collecting results
  # The goal is to ensure consistent content of numerical outputs

  results <- run_scm_collect(p2, env, ctrl)
  #  saveRDS(results, "tests/testthat/test_data/run_collect_tidy_2spp.rds")
  ref <- readRDS(file.path(rprojroot::find_testthat_root_file(), "test_data/run_collect_tidy_2spp.rds"))
  
  # check columns,should contain auxillary variables
  expect_equal(ncol(results$species), 15)
  expect_contains(names(results$species), c("competition_effect", "net_mass_production_dt")) 
  expect_contains(names(results), names(ref)[-1])
  
  expect_equal(results$steps$time, ref$time)
  expect_equal(results$n_spp, ref$n_spp)
  expect_equal(results$offspring_production, ref$offspring_production)
  expect_equal(results$p, ref$p)
  expect_equal(results$env$light_availability, ref$env$light_availability |> dplyr::select(names(results$env$light_availability)))

  # Need to maniupulate species object for effective comnparison
  v1 <- results$species |> dplyr::arrange(species, time, node) |> dplyr::filter(node != 142)
  
  v2 <- ref$species |>  tidyr::drop_na() |> dplyr::select(names(v1)) |>
    dplyr::arrange(species, time, node) |> dplyr::slice(1:nrow(v1))  
  expect_equal(v1, v2)
})
