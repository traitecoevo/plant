
context("Extrinsic drivers")

test_that("ExtrinsicDrivers can be instantiated", {
  
  expect_silent(drv <- ExtrinsicDrivers())

  expect_contains(names(drv), c('clone', 'evaluate', 'evaluate_range', 'get_names', 'initialize', 'set_constant', 'set_extrapolate', 'set_variable'))
  
  expect_silent(drv$set_constant("test", 1000))
 
  # extract
  expect_equal(drv$evaluate("test", 1), 1000)
  expect_equal(drv$evaluate("test", 10), 1000)

  # set new value
  expect_silent(drv$set_constant("test", 1.5))
  expect_equal(drv$evaluate("test", 1), 1.5)
  expect_equal(drv$evaluate("test", 10), 1.5)

  # set a bunch of values
  for(i in 1:10)
    expect_silent(drv$set_constant(letters[i], i))
  for (i in 1:10)
    expect_equal(drv$evaluate(letters[i], 1), i)

  # set new value
  for(i in 1:10)
    expect_silent(drv$set_constant(letters[i], i-1))
  for (i in 1:10)
    expect_equal(drv$evaluate(letters[i], 1), i-1)
 
  # expected erors
  expect_error(drv$evaluate("wrong", 1))
  
  # a function
  x <- 1:10
  y <- sin(x)
  expect_silent(drv$set_variable("sin", x,y))
  for (i in 1:10)
    expect_equal(drv$evaluate("sin", i), y[i])
})

test_that("ExtrinsicDrivers query methods", {

  drv <- ExtrinsicDrivers()
  drv$set_constant("a", 1)
  drv$set_variable("line", seq(0, 10), seq(0, 10))

  # get_names returns every active driver (order is unspecified)
  expect_setequal(drv$get_names(), c("a", "line"))

  # evaluate_range returns one value per query point, matching evaluate()
  pts <- c(0, 2.5, 7, 10)
  expect_equal(drv$evaluate_range("a", pts), rep(1, length(pts)))
  expect_equal(drv$evaluate_range("line", pts), pts)
  expect_equal(drv$evaluate_range("line", pts),
               vapply(pts, drv$evaluate, numeric(1), driver_name = "line"))

  # a linear spline interpolates between control points, not just at them
  expect_equal(drv$evaluate("line", 3.3), 3.3)
})

test_that("ExtrinsicDrivers extrapolation can be toggled", {

  drv <- ExtrinsicDrivers()
  drv$set_variable("line", seq(0, 10), seq(0, 10))

  # by default the spline does not extrapolate beyond its control points
  expect_error(drv$evaluate("line", 20))

  # once enabled, evaluation outside the range extends the trend
  expect_silent(drv$set_extrapolate("line", TRUE))
  expect_equal(drv$evaluate("line", 20), 20)
})

test_that("Variable birth rate is driven through the SCM", {
  # Closes the "scm untested" TODO on the variable extrinsic-driver path:
  # set a time-varying birth rate and confirm it flows through a full run.
  p0 <- scm_base_parameters("FF16")
  lma <- trait_matrix(0.0825, "lma")

  x_pts <- seq(0, 200)
  birth_rates <- list(sp1 = list(x = x_pts, y = 1 + 0.5 * sin(x_pts / 10)))

  p1 <- expand_parameters(lma, p0, FF16_hyperpar,
                          keep_existing_strategies = FALSE,
                          birth_rate_list = birth_rates)

  expect_true(p1$strategies[[1]]$is_variable_birth_rate)

  out <- run_scm(p1)
  drv <- out$patch$species[[1]]$extrinsic_drivers

  # the spline is active in the run: birth rate varies over time and tracks
  # the control points used to build it.
  expect_equal(drv$evaluate("birth_rate", 0), 1)
  expect_equal(drv$evaluate("birth_rate", 50), 1 + 0.5 * sin(50 / 10))
  expect_false(isTRUE(all.equal(drv$evaluate("birth_rate", 0),
                                drv$evaluate("birth_rate", 50))))

  # and the run itself produces a finite reproductive output
  expect_true(is.finite(out$net_reproduction_ratios))
})
