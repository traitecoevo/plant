

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

test_that("generate_strategy() gives a clear error when birth_rate$x is mis-specified (issue #381)", {
  p0 <- testing_parameters("FF16")
  lma <- trait_matrix(0.0825, "lma")
  max_life <- p0$max_patch_lifetime

  # x not starting at 0 → clear error
  expect_error(
    add_strategies(p0, lma, hyperpar = default_hyperpar("FF16"),
                   birth_rate = list(list(x = seq(1, max_life + 10), y = rep(1, max_life + 10))),
                   keep_existing = FALSE),
    regexp = "start at 0"
  )

  # x not extending past max_patch_lifetime → clear error
  expect_error(
    add_strategies(p0, lma, hyperpar = default_hyperpar("FF16"),
                   birth_rate = list(list(x = seq(0, max_life - 1), y = rep(1, max_life))),
                   keep_existing = FALSE),
    regexp = "max_patch_lifetime"
  )
})

test_that("Variable birth rate is driven through the SCM", {
  # Closes the "scm untested" TODO on the variable extrinsic-driver path:
  # set a time-varying birth rate and confirm it flows through a full run.
  p0 <- testing_parameters("FF16")
  lma <- trait_matrix(0.0825, "lma")

  x_pts <- seq(0, 200)
  birth_rates <- list(sp1 = list(x = x_pts, y = 1 + 0.5 * sin(x_pts / 10)))

  p1 <- add_strategies(p0, lma, hyperpar = default_hyperpar("FF16"), birth_rate = birth_rates, keep_existing = FALSE)

  expect_true(p1$strategies[[1]]$is_variable_birth_rate)

  out <- run(p1)$solver
  drv <- out$patch$species[[1]]$extrinsic_drivers

  # the spline is active in the run: birth rate varies over time and tracks
  # the control points used to build it.
  expect_equal(drv$evaluate("birth_rate", 0), 1)
  expect_equal(drv$evaluate("birth_rate", 50), 1 + 0.5 * sin(50 / 10))
  expect_false(isTRUE(all.equal(drv$evaluate("birth_rate", 0),
                                drv$evaluate("birth_rate", 50))))

  # and the run itself produces a finite reproductive output
  expect_true(is.finite(out$patch$net_reproduction_ratios))
})

# --- Time-varying climate drivers in TF24 ------------------------------------
#
# A step in incoming shortwave from 900 to 300 W m^-2 just after t = 1. The
# monotone spline gives a flat segment exactly between equal control points, so
# the stepped driver is identical to a constant 900 on [0, 1] and to a constant
# 300 from t = 1.1. The control points run well past the last time asked for,
# because an adaptive step may evaluate beyond it and a spline does not
# extrapolate by default.
swdown_step <- list(x = c(0, 0.5, 1, 1.1, 1.5, 2, 3, 5),
                    y = c(900, 900, 900, 300, 300, 300, 300, 300))

tf24_env_with_swdown <- function(swdown) {
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.2, 5))
  if (is.list(swdown)) {
    env$extrinsic_drivers_set_variable("SWdown", swdown$x, swdown$y)
  } else {
    env$extrinsic_drivers_set_constant("SWdown", swdown)
  }
  env
}

test_that("The stepped SWdown driver is flat either side of its step", {
  env <- tf24_env_with_swdown(swdown_step)
  expect_equal(env$extrinsic_drivers_evaluate_range("SWdown", seq(0, 1, by = 0.05)),
               rep(900, 21))
  expect_equal(env$extrinsic_drivers_evaluate_range("SWdown", seq(1.1, 5, by = 0.05)),
               rep(300, 79))
})

test_that("A lone TF24 individual follows a time-varying SWdown", {
  # A lone individual reads the atmospheric drivers at the runner's time. (Its
  # soil water is not integrated, so a rainfall driver would have no effect
  # here; that is a patch-level process.)
  grow <- function(swdown) {
    env <- tf24_env_with_swdown(swdown)
    env$set_fixed_environment(1.0, 150)
    ind <- Individual("TF24", "TF24_Env")(TF24_Model())
    res <- run(ind, times = c(0.5, 1, 2), env = env)
    res$state[, "height"]
  }
  h_hi <- grow(900)
  h_lo <- grow(300)
  h_step <- grow(swdown_step)

  # Less light, slower growth: the comparison means something.
  expect_true(all(h_lo < h_hi))

  # Before the step the stepped run sees the same driver as the constant 900
  # run. The solver's step grid near t = 1 can differ between the two, so this
  # is equal to within tolerance rather than identical.
  expect_equal(h_step[1:2], h_hi[1:2])

  # After it, growth slows: the stepped plant ends between the two constants.
  expect_lt(h_step[3], h_hi[3])
  expect_gt(h_step[3], h_lo[3])
})

test_that("A TF24 patch follows a time-varying SWdown", {
  # Height of the tallest node at each recorded time. Offspring production is
  # no use here: nothing is mature after two years, so it is below 1e-14 in
  # every run.
  tallest <- function(swdown) {
    p <- testing_parameters("TF24")
    p$max_patch_lifetime <- 2
    p <- add_strategies(p, trait_matrix(0.07, "lma"))
    d <- tidy(run(p, env = tf24_env_with_swdown(swdown), record = "steps"))$species
    tapply(d$height, d$time, max, na.rm = TRUE)
  }
  h_hi <- tallest(900)
  h_lo <- tallest(300)
  h_step <- tallest(swdown_step)

  # The node schedule fixes the recorded times, so the runs line up.
  expect_identical(names(h_step), names(h_hi))
  t <- as.numeric(names(h_hi))
  n <- length(t)

  expect_lt(h_lo[[n]], h_hi[[n]])

  # Before the step the stepped patch sees the same driver as the constant 900
  # patch.
  expect_equal(h_step[t <= 1], h_hi[t <= 1])

  # After it, growth slows: the stepped patch ends between the two constants.
  expect_lt(h_step[[n]], h_hi[[n]])
  expect_gt(h_step[[n]], h_lo[[n]])
})
