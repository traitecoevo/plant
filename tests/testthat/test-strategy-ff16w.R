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
  expect_equal(e$soil$states, 0)

  e$compute_rates(c(0.4))
  # 1 - 0.4 - 0.0*0.1 = 0.6
  expect_equal(e$soil$rates, 0.6) # default rainfall is now y = 1

  # Make it rain
  rain = list(
    x = seq(0, 9, 1),
    y = rep(5, 10)
  )
  e <- make_environment("FF16w", rainfall=rain)
  
  expect_equal(e$soil$states, 0)
  e$compute_rates(c(0.4))
  # 5 - 0.4 - 0.0*0.1
  expect_equal(e$soil$rates, 4.6)

  # Water logged
  e$set_soil_water_state(100)
  expect_equal(e$soil$states, 100)

  # Check construction
  e <- make_environment("FF16w",
                        soil_number_of_depths = 1,
                        rainfall = 10)

  expect_equal(e$soil$states, 0)
  e$compute_rates(c(0.4))
  # 10 - 0.4 - 0.0*0.1
  # weird because setting the soil water state to 100 doesn't affect
  # the rate at all since its still soil_layer 0
  expect_equal(e$soil$rates, 9.6)

  e <- make_environment("FF16w", soil_number_of_depths = 2)

  expect_equal(e$ode_size, 2)
  expect_equal(e$soil_number_of_depths, 2)
  expect_equal(e$soil$state_size, 2)

  expect_error(e <- make_environment("FF16w", soil_number_of_depths = 0),
               "FF16w Environment must have at least one soil layer")

  expect_error(e <- make_environment("FF16w",
                                     soil_number_of_depths = 1,
                                     soil_initial_state = c(1, 1)),
               "Not enough starting points for all layers")


})

test_that("Rainfall spline basic run", {
  # one species
  p0 <- scm_base_parameters("FF16w")
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16w_hyperpar,
                          mutant = FALSE, birth_rate_list = list(20))


  # init rainfall spline for env
  x <- seq(0, 110, 0.1)
  rain = list(
    x = x,
    y = 1000 + sin(x)
  )
  
  env <- make_environment("FF16w",
                          soil_number_of_depths = 10,
                          soil_initial_state = rep(1, 10),
                          rainfall = rain)

  ctrl <- scm_base_control()

  out <- run_scm(p1, env, ctrl)
  
  expect_equal(out$patch$environment$ode_size, 10)

  # This test only validates reproducibility across operating systems,
  # not any kind of ecological process. It should be replaced once the
  # water model is completed
  expect_equal(out$patch$environment$soil$rates,
               c(-0.9529, 0.2508, 1.4740, 5.1409, 13.5231, 28.4528, 49.8749,
                 74.9231, 98.4965, 115.2373), tolerance = 1e-5)

  # check the states are correct
  # again these values not yet verified
  expect_equal(out$patch$environment$soil$states,
               c(9939.809, 9877.551, 9803.061, 9691.901, 9496.920, 9152.641,
                 8594.141, 7785.160, 6740.445, 5528.322), tolerance = 1e-5)

  expect_equal(out$offspring_production, 16.88961056463, tolerance=5e-4)
  expect_equal(out$ode_times[c(10, 100)], c(0.000070, 4.101857), tolerance=1e-7)
})

test_that("Rainfall in collected output", {
  # one species
  p0 <- scm_base_parameters("FF16w")
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16w_hyperpar,
                          mutant = FALSE, birth_rate_list = list(20))
  
  
  env <- make_environment("FF16w",
                          soil_number_of_depths = 10,
                          soil_initial_state = rep(1, 10),
                          rainfall = 45)
  
  ctrl <- scm_base_control()
  
  collected <- run_scm_collect(p1, env, ctrl)
  
  expect_equal(collected$env[[1]]$rainfall, 45)
  expect_equal(collected$env[[142]]$rainfall, 45)
})
