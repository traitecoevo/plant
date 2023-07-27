context("fitness support")

test_that("positive_1d", {
  f <- function(x) -x^2 + 1
  tol <- 1e-8
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_lt(r[[1]], -1 + tol)
  expect_gt(r[[2]], 1 - tol)

  tol <- 1e-3
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_lt(r[[1]], -1 + tol)
  expect_gt(r[[2]], 1 - tol)

  expect_error(positive_1d(f, -2, 0.1, tol=tol), "no positive values")
})

test_that("bounds", {

  # First test object with infinite bounds
  expect_silent({ 
    bounds0 <- bounds_infinite("lma")
  })
  expect_is(bounds0, "matrix")
  expect_equal(bounds0,
    matrix(c(-Inf, Inf), nrow=1, dimnames =list("lma",c("lower", "upper")))
  )
  
  # Manually created bounds

  expect_silent(
    bounds_1d <- bounds(lma=c(0.01, 10))
  )
  expect_silent(
    bounds_2d <- bounds(lma=c(0.01, 10), rho=c(1, 1000))
  )

  expect_is(bounds_1d, "matrix")
  expect_equal(
    bounds_1d,
    matrix(c(0.01, 10), nrow = 1, dimnames = list("lma", c("lower", "upper")))
  )

  expect_is(bounds_2d, "matrix")
  expect_equal(
    bounds_2d,
    matrix(c(0.01, 10, 1, 1000), byrow=TRUE, nrow = 2, dimnames = list(c("lma", "rho"), c("lower", "upper")))
  )

  expect_silent(
    check_bounds(bounds_1d)
  )
  expect_silent(
    check_bounds(bounds_2d)
  )
  
  expect_error(
    check_bounds(c(0,1))
  )
  expect_error(
    check_bounds(matrix(c(0.01, 10, 1, 1000)))
  )
  expect_error(
    check_bounds(bounds0, finite=TRUE)
  )
  expect_silent(
    check_bounds(bounds0)
  )
  

  # Points lie within bounds
  expect_silent(
    check_point(0.02, bounds_1d)
  )
  expect_error(
    check_point(0.001, bounds_1d)
  )
})

test_that("fitness", {
  
  # Now solve actual bounds for viable fitness
  params <- scm_base_parameters("FF16")

  expect_silent({ 
    bounds1 <- viable_fitness(bounds_infinite("lma"), params)
  })

  expect_is(bounds1, "matrix")
  expect_equal(
      bounds1,
      matrix(c(0.02533822, 4.989169), nrow = 1, dimnames = list("lma", c("lower", "upper")))
  )
  expect_silent({
     bounds1b <- viable_fitness(bounds1, params)
  })
  expect_equal(
    bounds1b,
    matrix(c(0.02533822, 4.989169), nrow = 1, dimnames = list("lma", c("lower", "upper"))), 
    tolerance = 1e-4
  )
  # Now check edge cases

  ## Incorrect inputs
  expect_error(
    bounds2 <- viable_fitness(NA, params)
  )
  expect_error(
    bounds2 <- viable_fitness(c(-Inf, Inf), params)
  )
  ## starting point outside allowable area
  expect_error(
    bounds2 <- viable_fitness(bounds, x = -0.01, params)
  )
  ## 2D case currently failing
  expect_error(
    bounds2 <- viable_fitness(rbind(bounds0, bounds0, bounds0), params)
  )

  # max growth rate function
  expect_silent(
    fitness2 <- fundamental_fitness(trait_matrix(0.05, "lma"), params)
  )
  expect_equal(
    fitness2,
    8.372336,
    tolerance = 1e-4
  )

  # max fitness
  expect_silent(
    max_f <- max_fitness(bounds1, params, log_scale = TRUE)
  )

  expect_equal(max_f[1], 0.2985623, tolerance = 1e-4)
  expect_equal(attr(max_f, "fitness"), 11.09135, tolerance = 1e-4)

  # fitness landscape
  lma <- trait_matrix(seq_log_range(bounds1, 5), "lma")

  expect_silent({ 
    fitness <- fitness_landscape(lma, params)
  })

  comparison <- c(-0.000510755977418667, 10.3682397084861, 11.0785948241787, 9.84121044132885, -0.000147021768340749)

  expect_equal(
    fitness, comparison, 
    tolerance = 1e-4
  )
})

test_that("viable strategies", {
  params <- scm_base_parameters("FF16")
  params$max_patch_lifetime <- 60

  patch <- expand_parameters(trait_matrix(c(0.005, 0.03, 0.1), "lma"), params, birth_rate_list = c(1, 0.1, 0.001))

  ctrl = scm_base_control()

  expect_silent(
    ret <- check_inviable(patch, ctrl)
  )

  comparison <- c(0.0, 0.001386715, 3.150417828)
  expect_equal(as.numeric(ret), comparison, tolerance = 1e-4)
  expect_equal(attr(ret, "drop"), c(TRUE, FALSE, FALSE))

  # Change the extinction threshold
  ctrl$equilibrium_extinct_birth_rate <- 1
  expect_silent(
    ret <- check_inviable(patch, ctrl)
  )
  expect_equal(as.numeric(ret), comparison, tolerance = 1e-4)
  expect_equal(attr(ret, "drop"), c(TRUE, TRUE, FALSE))

})

