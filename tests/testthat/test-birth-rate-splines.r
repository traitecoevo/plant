context("Extrinisic drivers")

test_that("Can set birth rate splines correctly", {
  p0 <- scm_base_parameters("FF16")
  
  # two species
  lmas <- trait_matrix(c(0.0825, 0.125), "lma")
  
  # constant birth rates
  x <- seq(1, 200, len = 200)
  
  birth_rates <- list(
    species1 = list(x = x, y = 1 + sin(x)),
    species2 = list(x = x, y = 1 + cos(x))
  )
  
  
  p1 <- expand_parameters(lmas, p0, FF16_hyperpar, FALSE, birth_rates)
  
  # no longer stored in parameters
  p1$birth_rate
  p1$strategies[[1]]$birth_rate_x
  p1$strategies[[1]]$birth_rate_y
  p1$strategies[[1]]$is_variable_birth_rate
  
  
  # strategy
  s <- FF16_Species(p1$strategies[[1]]) # unable to initialise Interpolator
  
  s$extrinsic_drivers$evaluate("birth_rate", 200) # might need changes to RcppR6 API
  s$extrinsic_drivers()$evaluate_range("birth_rate", c(1, 2, 3))
  s$extrinsic_drivers()$get_names
  
  # trial run
  env <- make_environment("FF16")
  ctrl <- scm_base_control()
  
  out <- run_scm(p1, env, ctrl)

  # expect_something
})