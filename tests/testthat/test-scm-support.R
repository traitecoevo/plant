context("SCM support")


test_that("collect", {
  
  env <- Environment("FF16")
  ctrl <- scm_base_control()
  p0 <- scm_base_parameters("FF16")
  p0$disturbance_mean_interval <- 30.0
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, birth_rate_list = 1.0)

  expect_silent(res <- run_scm(p1, env, ctrl))

  expect_contains(
    names(res), cc("clone", "collect", "competition_effect_error", "complete", "history", "initialize", "net_reproduction_ratio_errors", "net_reproduction_ratio_for_species", "net_reproduction_ratios", "node_schedule", "ode_times", "offspring_production", "parameters", "patch", "reset", "run", "run_mutant", "run_next", "set_node_schedule_times", "time", "use_ode_times")
  )

  cmp_patch_density <- Weibull_Disturbance_Regime(p1$max_patch_lifetime)$density(res$time)
  expect_equal(res$patch_density, cmp_patch_density)

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

  results <- run_scm_collect(p2, env, ctrl)
  
  expect_equal(ncol(results$species), 14)
  
  v <- c("competition_effect", "net_mass_production_dt")
  expect_contains(names(results$species), c("competition_effect", "net_mass_production_dt")) 
#  saveRDS(results, "tests/testthat/test_data/run_collect_tidy_2spp.rds")
  
  ref <- readRDS(file.path(rprojroot::find_testthat_root_file(), "test_data/run_collect_tidy_2spp.rds"))
  expect_contains(names(results), names(ref))
  
  expect_equal(results$time, ref$time)
  expect_equal(results$n_spp, ref$n_spp)
  expect_equal(results$offspring_production, ref$offspring_production)
  expect_equal(results$p, ref$p)
  expect_equal(results$env$light_availability, ref$env$light_availability |> dplyr::select(-patch_density))

  v1 <- results$species |> dplyr::arrange(species, time, node) 
  
  v2 <- ref$species |>
    tidyr::drop_na() |>
    dplyr::select(names(v1)) |>
    dplyr::arrange(species, time, node) |> dplyr::slice(1:nrow(v2))
  
  # failing here
  expect_equal(v1, v2)
  

})
