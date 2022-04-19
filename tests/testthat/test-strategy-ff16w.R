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
  
  e$compute_rates()
  expect_equal(e$soil$rates, 1) # default rainfall is now y = 1
  
  # Make it rain
  x <- seq(0, 9, 1)
  y <- rep(5, 10)
  e$set_extrinsic_driver("rainfall", x, y)
  expect_equal(e$soil$states, 0)
  e$compute_rates()
  expect_equal(e$soil$rates, 5)
  
  # Water logged
  e$set_soil_water_state(100)
  expect_equal(e$soil$states, 100)
  
  # Check construction
  e <- make_environment("FF16w", 
                        soil_number_of_depths = 1,
                        rainfall = 10)
  
  expect_equal(e$soil$states, 0)
  e$compute_rates()
  expect_equal(e$soil$rates, 10)

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
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16w_hyperpar,FALSE,20)
  #p1$birth_rate <- 20
  
  env <- make_environment("FF16w", 
                          soil_number_of_depths = 10,
                          soil_initial_state = rep(1, 10))
  
  # init rainfall spline for env
  x <- seq(0, 110, 0.1)
  integrand <- function(x) {x^2}
  y <- integrand(x)
  env$set_extrinsic_driver("rainfall", x, y)
  
  ctrl <- scm_base_control()
  
  out <- run_scm(p1, env, ctrl)
  
  expect_equal(out$patch$environment$ode_size, 10)
  
  # model is currently 1/z * int(0, t)(f(x)), where z is depth, 
  # x is time, f is spline, t is number of years
  
  # check the rates are correct, ie 105^2/depth
  expect_equal(out$patch$environment$soil$rates, 
               sapply(seq(1, 10), function(x) { (out$time^2)/x }),
               tolerance = 1e-7)
  
  # check the states are correct
  # model is currently 1/z * int(0, t)(f(x)), where z is depth, 
  # x is time, f is spline, t is number of years
  expect_equal(out$patch$environment$soil$states,
               sapply(seq(1, 10), function(i) { (1/i) * integrate(integrand, 0, out$time)$value }),
               tolerance = 1e-5)
  
  expect_equal(out$offspring_production, 16.88946, tolerance=1e-5)
  expect_equal(out$ode_times[c(10, 100)], c(0.000070, 4.216055), tolerance=1e-7)
})