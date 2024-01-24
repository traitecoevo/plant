# Built from  tests/testthat/test-strategy-ff16.R on Mon Jul 19 11:01:04 2021 using the scaffolder, from the strategy:  FF16
context("Strategy-FF16w")

test_that("FF16w Environment", {
  control <- Control()
  e <- make_environment("FF16w")

  # One layer by default
  expect_equal(e$ode_size, 1)
  expect_equal(e$soil_number_of_depths, 1)
  expect_equal(e$soil$state_size, 1)

  # Initialised with no water, no inflow
  e$set_soil_water_state(0)
  expect_equal(e$soil$states,0)
  
  e$extrinsic_drivers$evaluate("rainfall", 1)
  
  e$compute_rates(c(0.4))
  # 1 - 0.4  = 0.6
  expect_equal(e$soil$rates, 0.6) # default rainfall is now y = 1

  # Make it rain
  rain = list(
    x = seq(0, 9, 1),
    y = rep(0.6, 10)
  )
  e <- make_environment("FF16w", rainfall=rain, soil_initial_state = 0)

  expect_equal(e$soil$states, 0)
  e$compute_rates(c(0.3))
  # 0.6 - 0.3 
  expect_equal(e$soil$rates, 0.3)

  # Water logged
  e$set_soil_water_state(100)
  expect_equal(e$soil$states, 100)

  # Check construction
  e <- make_environment("FF16w",
                        soil_number_of_depths = 1,
                        rainfall = 10,
                        soil_initial_state = 0)

  expect_equal(e$soil$states, 0)
  
  e$compute_rates(c(0.4))
  # 10 - 0.4 - 0.0*0.1
  # weird because setting the soil water state to 100 doesn't affect
  # the rate at all since its still soil_layer 0
  expect_equal(e$soil$rates, 9.6)

  e <- make_environment("FF16w", soil_number_of_depths = 2, soil_initial_state = c(0,0))

  expect_equal(e$ode_size, 2)
  expect_equal(e$soil_number_of_depths, 2)
  expect_equal(e$soil$state_size, 2)

  expect_error(e <- make_environment("FF16w", soil_number_of_depths = 0),
               "FF16w Environment must have at least one soil layer")

  expect_error(e <- make_environment("FF16w",
                                     soil_number_of_depths = 1,
                                     soil_initial_state = c(1, 1)),
               "Not enough starting points for all layers")

  
  # Check that having a slightly full soil layer reduces the rate of recharge
  e <- make_environment("FF16w",
                        soil_number_of_depths = 1,
                        rainfall = 1,
                        soil_initial_state = 0.4)
  
  e$compute_rates(c(0.4))
  saturated_rate <- e$soil$rates
  
  e$set_soil_water_state(0)
  e$compute_rates(c(0.4))
  unsaturated_rate <- e$soil$rates
  
  expect_true(saturated_rate < unsaturated_rate)
  
  e <- FF16w_make_environment(ca = 42,
                        atm_vpd = 1.5,
                        leaf_temp = 22.5,
                        atm_o2_kpa = 123,
                        atm_kpa = 2.1)
  
  expect_equal(e$get_atm_kpa(), 2.1)
  expect_equal(e$get_atm_o2_kpa(), 123)
  expect_equal(e$get_leaf_temp(), 22.5)
  expect_equal(e$get_atm_vpd(), 1.5)
  expect_equal(e$get_ca(), 42)
  
  atm_vpd = list(
    x = seq(0, 9, 1),
    y = seq(0.6, 1, length.out = 10)
  )
  
  e <- FF16w_make_environment(atm_vpd = atm_vpd)
  e$extrinsic_drivers$evaluate("atm_vpd", 2)
  
  expect_equal(e$extrinsic_drivers$evaluate("atm_vpd", 2), 0.6888889)

})

test_that("Rainfall spline basic run", {
  # one species
  p0 <- scm_base_parameters("FF16w")
  p0$max_patch_lifetime <- 10
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16w_hyperpar,
                           birth_rate_list = list(10))


  # init rainfall spline for env
  x <- seq(0, 10, length.out = 1000)
  rain = list(
    x = x,
    y = 2 + sin(x)
  )

  env <- make_environment("FF16w",
                          soil_number_of_depths = 1,
                          soil_initial_state = rep(0.4),
                          rainfall = rain)

  ctrl <- scm_base_control()
  ctrl$ci_niter = 1000
  ctrl$GSS_tol_abs <- 1e-1
  ctrl$ci_abs_tol <- 1e-8
  ctrl$vulnerability_curve_ncontrol <- 100

  out <- run_scm(p1, env, ctrl)

  expect_equal(out$patch$environment$ode_size, 1)

  # This test only validates reproducibility across operating systems,
  # not any kind of ecological process. It should be replaced once the
  # water model is completed
  expect_equal(out$patch$environment$soil$rates,
               c(-0.2374372), tolerance = 1e-5)

  # check the states are correct
  # again these values not yet verified
  expect_equal(out$patch$environment$soil$states,
               c(0.334688), tolerance = 1e-5)

  expect_equal(out$offspring_production, 0.0003052265, tolerance=5e-4)
  expect_equal(out$ode_times[c(1,5,10)], c(0, 2e-05, 7e-05), tolerance=1e-7)
})

test_that("Rainfall in collected output", {
  # one species
  p0 <- scm_base_parameters("FF16w")
  
  p0$max_patch_lifetime <- 7.01036

  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16w_hyperpar,
                           birth_rate_list = list(100))
  
  
  env <- make_environment("FF16w",
                          soil_number_of_depths = 1,
                          soil_initial_state = rep(0.4, 1),
                          rainfall = 1)
  
  ctrl <- scm_base_control()
  ctrl$ci_niter = 1000
  ctrl$GSS_tol_abs <- 1e-1
  ctrl$ci_abs_tol <- 1e-8
  ctrl$vulnerability_curve_ncontrol <- 100
  
  system.time(collected <- run_scm_collect(p1, env, ctrl))
  
  expect_equal(collected$env[[1]]$rainfall, 45)
  expect_equal(collected$env[[142]]$rainfall, 45)
})
