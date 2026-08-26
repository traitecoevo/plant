
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()


test_that("Environment-TF24 drivers", {
  
  
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
  expect_equal(env$extrinsic_drivers_evaluate("atm_kpa", 0), 101.3)
  expect_equal(env$extrinsic_drivers_evaluate("atm_kpa", 10), 101.3)

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
  
  theta_sat <- 0.428
  a_psi <- 1.78e3
  
  expect_silent(env$set_soil_water_state(rep(theta_sat,num_depths)))
  expect_equal(env$get_soil_water_state(), rep(theta_sat,num_depths))
  expect_equal(env$get_soil_water_state_cumulative_flux(), rep(0,5))  
  
  expect_equal(env$psi_from_soil_moist(theta_sat), a_psi/1e6)
  
  # test handy wrappers
  expect_equal(env$get_PPFD(), env$extrinsic_drivers_evaluate("PPFD", 1))
  expect_equal(env$get_atm_vpd(), env$extrinsic_drivers_evaluate("atm_vpd", 1))
  expect_equal(env$get_ca(), env$extrinsic_drivers_evaluate("ca", 1))
  expect_equal(env$get_leaf_temp(), env$extrinsic_drivers_evaluate("leaf_temp", 1))
  expect_equal(env$get_atm_o2_kpa(), env$extrinsic_drivers_evaluate("atm_o2_kpa", 1))
  expect_equal(env$get_atm_kpa(), env$extrinsic_drivers_evaluate("atm_kpa", 1))

  
})

test_that("Environment-TF24 soil layers", {

  env <- Environment("TF24")
  theta_sat <- 0.428
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

test_that("Environment-TF24 soil moisture and potential invert each other", {

  env <- Environment("TF24")
  theta_sat <- 0.428
  theta_residual <- 1e-2

  theta <- seq(theta_residual, theta_sat, length.out = 200)[-1]
  psi <- vapply(theta, env$psi_from_soil_moist, numeric(1))

  # psi_from_soil_moist caps its output, so the inverse can only recover the
  # moistures whose potential sits below that ceiling.
  below_ceiling <- psi < max(psi)
  expect_true(sum(below_ceiling) > 100)

  theta_back <- vapply(psi[below_ceiling], env$soil_moist_from_psi, numeric(1))
  expect_equal(theta_back, theta[below_ceiling], tolerance = 1e-12)
  expect_true(all(theta_back <= theta_sat))
})

test_that("Environment-TF24 running soil moisture profile", {

  env <- Environment("TF24")
  # get list of extrinsic drivers for the environment

  # default values
  expect_equal(env$soil_moist_sat, 0.428)
  expect_equal(env$K_sat, 163.0411)
  expect_equal(env$a_psi, 1.78e3)
  expect_equal(env$n_psi, 6.57)
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
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"))
  out <- run_scm(p1, env)
  expect_equal(out$patch$environment$n_psi, 4)
  expect_equal(out$patch$environment$b_infil, 5)
  
  env <- Environment("TF24")
  env$set_soil_number_of_depths(1)
  env$set_soil_water_state(0.1)
  p0 <- scm_base_parameters("TF24")
  out <- run_scm(p1, env, collect = TRUE)
  expect_equal(length(unique(out$env$soil_depth$soil_depth)), 1)

  
  # check conservation of water for 1 layer
  
  depth <- out$env$soil_depth$soil_depth[1]
  out$env$soil_moist_cumulative_flux %>%
    dplyr::mutate(sum_runoff = sum_rainfall - sum_infiltration) -> cumulative_fluxes
  
  out$env$soil_moist %>%
    dplyr::mutate(soil_moist_mm = soil_moist*depth) -> water_storage
  
  cumulative_fluxes %>%
    dplyr::left_join(water_storage,
                     by = c("time", "step", "patch_density")) %>%
    dplyr::slice(c(1,nrow(.))) %>%
    dplyr::mutate(init_soil_moist_mm = soil_moist_mm[1]) %>%
    dplyr::slice_tail(n = 1) %>%
    # water from total rainfdall over period and initial storage
    dplyr::mutate(total_moisture_start = sum_rainfall + init_soil_moist_mm) %>%
    # water in storage at end plus water lost to bottom drainage and runoff
    dplyr::mutate(total_moisture_end = soil_moist_mm + sum_drainage + sum_runoff) -> water_conservation
  
  expect_equal(water_conservation$total_moisture_start, water_conservation$total_moisture_end)
  
})

test_that("Environment-TF24 allows per-layer soil parameters", {
  env <- TF24_Environment()

  expect_silent(env$set_soil_parameters(3, NULL, NULL, NULL, NULL))
  expect_equal(env$get_soil_number_of_depths(), 3)
  expect_silent(env$set_soil_water_state(rep(0.2, 3)))

  # Provided parameter vectors must match number of soil layers.
  expect_error(
    env$set_soil_parameters(3, c(0.4, 0.5), c(10, 20, 30), c(1000, 2000, 3000), c(1, 1, 1)),
    "soil_moist_sat"
  )

  expect_silent(
    env$set_soil_parameters(3, c(0.4, 0.5, 0.6), c(10, 20, 30), c(1000, 2000, 3000), c(1, 1, 1))
  )
  expect_silent(env$set_soil_water_state(rep(0.2, 3)))

  # The per-layer parameterised environment should run rate calculations
  # without errors for the configured number of soil layers.
  expect_silent(env$compute_rates(rep(0, 3)))

  # If parameter vectors are not provided (NULL), defaults are replicated
  # from the scalar environment parameters across layers.
  env2 <- TF24_Environment()
  env2$soil_moist_sat <- 0.5
  env2$a_psi <- 1e3
  env2$n_psi <- 1
  expect_silent(env2$set_soil_parameters(2, NULL, NULL, NULL, NULL))
  expect_silent(env2$set_soil_water_state(rep(0.2, 2)))
  expect_silent(env2$compute_rates(rep(0, 2)))
  expect_equal(length(env2$get_soil_water_state()), 2)
})

test_that("Environment-TF24 clear returns the soil to its starting state", {
  env <- Environment("TF24")
  n <- env$get_soil_number_of_depths()

  expect_silent(env$clear())
  expect_equal(env$get_soil_water_state(), rep(0.428 * 0.5, n))
  expect_equal(env$get_soil_water_state_cumulative_flux(), rep(0, 5))

  # clear() returns the state last set, not the constructed default.
  start <- seq(0.30, by = 0.01, length.out = n)
  env$set_soil_water_state(start)
  env$clear()
  expect_identical(env$get_soil_water_state(), start)
  expect_equal(env$get_soil_water_state_cumulative_flux(), rep(0, 5))

  # After a changed layer count, the starting state is the one set since.
  env$set_soil_parameters(3, NULL, NULL, NULL, NULL)
  env$set_soil_water_state(c(0.2, 0.25, 0.3))
  env$clear()
  expect_identical(env$get_soil_water_state(), c(0.2, 0.25, 0.3))
})

test_that("Environment-TF24 set_soil_parameters validates each parameter length", {
  env <- TF24_Environment()

  expect_error(
    env$set_soil_parameters(3, c(0.4, 0.5), NULL, NULL, NULL),
    "soil_moist_sat"
  )
  expect_error(
    env$set_soil_parameters(3, NULL, c(10, 20), NULL, NULL),
    "K_sat"
  )
  expect_error(
    env$set_soil_parameters(3, NULL, NULL, c(1000, 2000), NULL),
    "a_psi"
  )
  expect_error(
    env$set_soil_parameters(3, NULL, NULL, NULL, c(1, 1)),
    "n_psi"
  )
})
