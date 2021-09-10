# Built from  tests/testthat/test-strategy-ff16.R on Mon Jul 19 11:01:04 2021 using the scaffolder, from the strategy:  FF16
context("Strategy-Water")

test_that("Water Environment", {
  control <- Control()
  e <- FF16_Environment(control)

  # Empty by default
  expect_equal(e$ode_size, 0)
  expect_equal(e$soil$state_size, 0)
  
  # Add one layer
  control$soil_number_of_depths <- 1
  e <- FF16_Environment(control)
  
  expect_equal(e$ode_size, 1)
  expect_equal(e$soil$state_size, 1)
  
  # Initialised with no water, no inflow
  expect_equal(e$soil$states, 0)
  
  e$compute_rates()
  expect_equal(e$soil$rates, 0)
  
  # Make it rain
  control$soil_infiltration_rate <- 10
  e <- FF16_Environment(control)
  
  expect_equal(e$soil$states, 0)
  e$compute_rates()
  expect_equal(e$soil$rates, 10)
})

test_that("Basic run", {
  # one species
  p0 <- scm_base_parameters("Water")
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, Water_hyperpar,FALSE)


  p1$birth_rate <- 20

  ctrl <- scm_base_control()
  ctrl$soil_infiltration_rate <- 1
  ctrl$soil_number_of_depths <- 10
  
  out <- run_scm(p1, ctrl)
  
  expect_equal(out$patch$environment$ode_size, 10)

  expect_equal(out$patch$environment$soil$rates, 
               c(1.00, 0.50, 0.33, 0.25, 0.20, 0.17, 0.14, 0.13, 0.11, 0.10),
               tolerance = .01)
  
  expect_equal(out$patch$environment$soil$states,
               c(105, 52, 35, 26, 21, 17, 15, 13, 11, 10),
               tolerance = 0.1)
  
  expect_equal(out$offspring_production, 16.88946, tolerance=1e-5)
  expect_equal(out$ode_times[c(10, 100)], c(0.000070, 4.216055), tolerance=1e-5)
})
