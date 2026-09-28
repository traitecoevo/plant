
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  e <- environment_types[[x]]


  test_that("Empty environment", {
    p <- Parameters(x, e)()
    env <- Environment(x)

    ## At this point, we should have full canopy openness, partly because
    ## the spline is just not constructed.
    expect_equal(env$get_environment_at_height(0), 1.0)
    expect_equal(env$get_environment_at_height(100), 1.0)

    ## The field answers queries before any build, from a flat interpolant over
    ## [0, 1]. Its knots are the first column of the state matrix.
    st <- env$light_availability$state
    expect_equal(nrow(st), 3)
    expect_equal(st[, 1], c(0, 0.5, 1))
    expect_equal(st[, 2], rep(1, 3))
    expect_equal(st[, 3], rep(0, 3))
  })

  test_that("Manually set environment", {
    env <- Environment(x)
    ## Now, set the light environment.
    hh <- seq(0, 10, length.out=101)
    light_env <- function(x) {
      exp(x/(max(hh)*2)) - 1 + (1 - (exp(.5) - 1))/2
    }
    ee <- light_env(hh)
    ## The field carries a slope per knot as well as a value, so it is set from
    ## heights, values and slopes rather than from a value-fitted interpolant.
    slope <- exp(hh/(max(hh)*2)) / (max(hh)*2)
    ## light_availability is a field, so reading it copies: mutate the copy and
    ## assign it back, or the environment keeps the field it already had.
    field <- env$light_availability
    field$init_interpolators(c(hh, ee, slope))
    env$light_availability <- field

    st <- env$light_availability$state
    expect_identical(st[, 1], hh)
    expect_identical(st[, 2], ee)

    ## A Hermite reproduces the value and the slope it was given at every knot,
    ## and interpolates between them to the accuracy of the data it was built on.
    expect_equal(sapply(hh, env$light_availability$get_value_at_height), ee)
    hmid <- (hh[-1] + hh[-length(hh)])/2
    expect_equal(sapply(hmid, env$light_availability$get_value_at_height),
                 light_env(hmid), tolerance = 1e-6)
  })
}

test_that("resource spline value is floored at zero (#253)", {
  # The cubic light spline can undershoot below zero between knots (notably the
  # K93 light availability at high k_I: a sharp Beer's-law drop to a low, flat
  # value). Build such a spline directly and confirm that
  # ResourceSpline::get_value_at_height never returns a negative resource
  # availability, while staying a no-op wherever the spline is non-negative.
  env <- Environment("K93")

  # A Hermite between non-negative knots undershoots only if its slopes are steep
  # enough to carry it below zero between them, so the slopes are what builds the
  # scenario here.
  hh <- c(0, 1, 2, 3, 4, 5)
  ee <- c(1, 1, 1, 0.02, 0.02, 0.02)
  slope <- c(0, 0, -12, -12, 0, 0)
  field <- env$light_availability
  field$init_interpolators(c(hh, ee, slope))
  env$light_availability <- field

  grid    <- seq(0, 5, by = 0.05)
  floored <- vapply(grid, env$light_availability$get_value_at_height, numeric(1))

  # The floor genuinely fires, otherwise the test proves nothing: the polynomial
  # goes negative between knots and the accessor returns exactly zero there.
  expect_true(any(floored == 0))
  # ... and it never returns a negative resource availability.
  expect_true(all(floored >= 0))
  # Away from the undershoot the knots are reproduced untouched.
  expect_equal(vapply(hh, env$light_availability$get_value_at_height, numeric(1)), ee)
})
