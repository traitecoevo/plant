context("Extrinisic drivers")

test_that("Can set birth rate splines correctly", {
  p0 <- scm_base_parameters("FF16")
  
  # two species
  lmas <- trait_matrix(c(0.0825, 0.125), "lma")
  
  # constant birth rates
  x <- seq(0, 200)
  
  birth_rates <- list(
    species1 = list(x = x, y = 1 + sin(x)),
    species2 = list(x = 1, y = 2) # probably a better way to set a constant birth_rate since we only use y
  )
  
  p1 <- expand_parameters(lmas, p0, FF16_hyperpar, FALSE, birth_rates)
  
  # no longer stored in parameters
  expect_null(p1$birth_rate)
  expect_equal(p1$strategies[[1]]$birth_rate_x, x)
  expect_equal(p1$strategies[[1]]$birth_rate_y,  1 + sin(x))
  expect_true(p1$strategies[[1]]$is_variable_birth_rate)
  expect_false(p1$strategies[[2]]$is_variable_birth_rate)
  
  
  # strategy 1 (variable)
  s1 <- FF16_Species(p1$strategies[[1]])

  k = c(1, 2, 3)
  expect_equal(s1$extrinsic_drivers()$evaluate("birth_rate", 200), 1 + sin(200))
  expect_equal(s1$extrinsic_drivers()$evaluate_range("birth_rate", k), 1 + sin(k))
  expect_equal(s1$extrinsic_drivers()$get_names(), c("birth_rate"))
  
  # strategy 2 (constant)
  s2 <- FF16_Species(p1$strategies[[2]])
  
  expect_equal(s2$extrinsic_drivers()$evaluate("birth_rate", 200), 2)
  expect_equal(s2$extrinsic_drivers()$evaluate_range("birth_rate", k), c(2, 2, 2))
  expect_equal(s2$extrinsic_drivers()$get_names(), c("birth_rate"))
  
  # trial run
  env <- make_environment("FF16")
  ctrl <- scm_base_control()
  
  out <- run_scm(p1, env, ctrl)

  # expect_something
})