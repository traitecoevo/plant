
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

context("Environment-TF24")

test_that("Environment-TF24 drivers", {
  
  context("TF24-Env-ExtrinsicDrivers")
  
  env <- Environment("TF24")
  # get list of extrinsic drivers for the environment

  expect_contains(env$extrinsic_drivers_get_names(), c("PPFD", "rainfall", "leaf_temp", "atm_o2_kpa", "atm_kpa", "ca", "atm_vpd"))
  
  # test default values - check at two values of second argument (should give same result)
  expect_equal(env$extrinsic_drivers_evaluate("PPFD", 0), 1800)
  expect_equal(env$extrinsic_drivers_evaluate("PPFD", 10), 1800)
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", 0), 1)
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", 10), 1)
  expect_equal(env$extrinsic_drivers_evaluate("atm_vpd", 0), 1)
  expect_equal(env$extrinsic_drivers_evaluate("atm_vpd", 10), 1)
  expect_equal(env$extrinsic_drivers_evaluate("ca", 0), 40)
  expect_equal(env$extrinsic_drivers_evaluate("ca", 10), 40)
  expect_equal(env$extrinsic_drivers_evaluate("PPFD", 0), 1800)
  expect_equal(env$extrinsic_drivers_evaluate("PPFD", 10), 1800)
  expect_equal(env$extrinsic_drivers_evaluate("atm_kpa", 0), 100.5)
  expect_equal(env$extrinsic_drivers_evaluate("atm_kpa", 10), 100.5)

  # test updating values
  v <- 200
  expect_silent(env$extrinsic_drivers_set_constant("rainfall", v))
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", 100), v)
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", 10000000), v)

  ## a function (simple quadratic)
  x <- seq(-10, 10, 0.41)
  y <- x^2
  PPFD <- 10
  
  env <- Environment("TF24")
  expect_silent(env$extrinsic_drivers_set_constant("PPFD", PPFD))
  expect_equal(env$extrinsic_drivers_evaluate("PPFD", 2), PPFD)
  
  # interpolated points
  expect_silent(env$extrinsic_drivers_set_variable("rainfall", x, y))
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", 2), 4)
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", -2), 4)
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", 3), 9, tolerance=1e-7)
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", -3), 9, tolerance=1e-7)
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", 5.5), 30.25, tolerance=1e-7)
  expect_equal(env$extrinsic_drivers_evaluate("rainfall", -5.5), 30.25, tolerance=1e-7)
  
  ## interpolated range of points
  expect_equal(env$extrinsic_drivers_evaluate_range("rainfall", c(-7, 1, 7.8345)), c(49, 1, 61.37939025), tolerance=1e-6)
  
  # test how many soil depths
  num_depths <- env$get_soil_number_of_depths()
  expect_equal(num_depths, 5)
  
  theta_sat <- 0.453
  a_psi <- 8.7
  
  expect_silent(env$set_soil_water_state(rep(theta_sat,num_depths)))
  expect_equal(env$get_soil_water_state(), rep(theta_sat,num_depths))
  expect_equal(env$get_soil_water_state_cumulative_flux(), rep(0,3))  
  
  expect_equal(env$psi_from_soil_moist(theta_sat), a_psi)
  
  # test handy wrappers
  expect_equal(env$get_PPFD(), env$extrinsic_drivers_evaluate("PPFD", 1))
  expect_equal(env$get_atm_vpd(), env$extrinsic_drivers_evaluate("atm_vpd", 1))
  expect_equal(env$get_ca(), env$extrinsic_drivers_evaluate("ca", 1))
  expect_equal(env$get_leaf_temp(), env$extrinsic_drivers_evaluate("leaf_temp", 1))
  expect_equal(env$get_atm_o2_kpa(), env$extrinsic_drivers_evaluate("atm_o2_kpa", 1))
  expect_equal(env$get_atm_kpa(), env$extrinsic_drivers_evaluate("atm_kpa", 1))

  
})

test_that("Environment-TF24 soil layers", {
  context("TF24-Env-Soil water")

  env <- Environment("TF24")
  theta_sat <- 0.453
  num_depths <- env$get_soil_number_of_depths()
  
  # get list of extrinsic drivers for the environment

  # default value is 1
  expect_equal(env$get_soil_number_of_depths(), 5)
  expect_equal(env$get_soil_water_state(), rep(theta_sat/2,num_depths))

  expect_silent(env$set_soil_water_state(rep(0.5,5)))
  expect_equal(env$get_soil_water_state(), rep(0.5,5))
  # should error when passed a vector that is too long
  expect_error(env$set_soil_water_state(c(0.5, 0.4)))

  # resize
  layers <- 10
  expect_silent(env$set_soil_number_of_depths(layers))
  expect_equal(env$get_soil_number_of_depths(), layers)
  expect_equal(env$get_soil_water_state(), rep(0, layers))

  expect_silent(env$set_soil_water_state(1:10))
  expect_equal(env$get_soil_water_state(), 1:10)
  # should error when passed a vector that is too long
  expect_error(env$set_soil_water_state(c(0.5, 0.4)))
 })

test_that("Environment-TF24 running soil moisture profile", {
  context("TF24-Env-parameters")

  env <- Environment("TF24")
  # get list of extrinsic drivers for the environment

  # default values
  expect_equal(env$soil_moist_sat, 0.453)
  expect_equal(env$K_sat, 440.628)
  expect_equal(env$a_psi, 8.7)
  expect_equal(env$n_psi, 4.8)
  expect_equal(env$a_infil, 1)
  expect_equal(env$b_infil, 8)

  # set values
  expect_silent(env$soil_moist_sat <- 1)
  expect_silent(env$K_sat <- 2)
  expect_silent(env$a_psi <- 3)
  expect_silent(env$n_psi <- 4)
  expect_silent(env$a_infil <- 0)
  expect_silent(env$b_infil <- 5)

  expect_equal(env$soil_moist_sat, 1)
  expect_equal(env$K_sat, 2)
  expect_equal(env$a_psi, 3)
  expect_equal(env$n_psi, 4)
  expect_equal(env$a_infil, 0)
  expect_equal(env$b_infil, 5)

  # check values from above are inherited when apssed into scm
  p0 <- scm_base_parameters("TF24")
  p0$max_patch_lifetime <- 0.01  
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0)
  out <- run_scm(p1, env)
  expect_equal(out$patch$environment$n_psi, 4)
  expect_equal(out$patch$environment$b_infil, 5)
  
  env <- Environment("TF24")
  env$set_soil_number_of_depths(1)
  env$set_soil_water_state(0.1)
  p0 <- scm_base_parameters("TF24")
  out <- run_scm_collect(p1, env)
  expect_equal(length(unique(out$env$soil_depth$soil_depth)), 1)

  
  # check conservation of water for 1 layer
  
  depth <- out$env$soil_depth$soil_depth[1]
  out$env$soil_moist_cumulative_flux %>%
    dplyr::mutate(sum_runoff = sum_rainfall - sum_infiltration) -> cumulative_fluxes
  
  out$env$soil_moist %>%
    dplyr::mutate(soil_moist_mm = soil_moist*depth) -> water_storage
  
  cumulative_fluxes %>%
    dplyr::left_join(water_storage) %>%
    dplyr::slice(c(1,nrow(.))) %>%
    dplyr::mutate(init_soil_moist_mm = soil_moist_mm[1]) %>%
    dplyr::slice_tail(n = 1) %>%
    # water from total rainfdall over period and initial storage
    dplyr::mutate(total_moisture_start = sum_rainfall + init_soil_moist_mm) %>%
    # water in storage at end plus water lost to bottom drainage and runoff
    dplyr::mutate(total_moisture_end = soil_moist_mm + sum_drainage + sum_runoff) -> water_conservation
  
  expect_equal(water_conservation$total_moisture_start, water_conservation$total_moisture_end)
  
})
