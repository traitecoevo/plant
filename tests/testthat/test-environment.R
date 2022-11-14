
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  e <- environment_types[[x]]

  context(sprintf("Environment-%s",x))

  test_that("Empty environment", {
    p <- Parameters(x, e)()
    env <- make_environment(x)

    ## At this point, we should have full canopy openness, partly because
    ## the spline is just not constructed.
    expect_equal(env$canopy_openness(0), 1.0)
    expect_equal(env$canopy_openness(100), 1.0)

    spline <- env$canopy$canopy_interpolator
    expect_equal(spline$size, 33)
    expect_equal(spline$x, seq(0,1, length.out=33))
  })

  test_that("Manually set environment", {
    env <- make_environment(x)
    ## Now, set the light environment.
    hh <- seq(0, 10, length.out=101)
    light_env <- function(x) {
      exp(x/(max(hh)*2)) - 1 + (1 - (exp(.5) - 1))/2
    }
    ee <- light_env(hh)
    interplator <- Interpolator()
    interplator$init(hh, ee)

    ## And set it
    env$canopy$canopy_interpolator <- interplator

    expect_identical(env$canopy$canopy_interpolator$xy, interplator$xy)

    hmid <- (hh[-1] + hh[-length(hh)])/2
    expect_identical(sapply(hmid, env$canopy$canopy_interpolator$eval), sapply(hmid, interplator$eval))
  })
}

test_that("FF16 rainfall spline", {
  
  context("Rainfall-FF16")
  
  env <- make_environment("FF16w")
  # get list of extrinsic drivers for the environment
  expect_equal(env$extrinsic_drivers$get_names(), c("rainfall"))
  
  # test extrapolation on default spline of y = 1
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", 100), 1)
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", 10000000), 1)
  
  # test extrapolation on spline of y = 5.613432
  env <- make_environment("FF16w", rainfall=5.613432)
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", 100), 5.613432)
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", 10000000), 5.613432)
  
  ## simple quadratic
  x <- seq(-10, 10, 0.41)
  quadratic_rain <- list(
    x = x,
    y = x^2
  )
  env <- make_environment("FF16w", rainfall=quadratic_rain) # overwrites previously created spline
  
  # interpolated points
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", 2), 4)
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", -2), 4)
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", 3), 9, tolerance=1e-7)
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", -3), 9, tolerance=1e-7)
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", 5.5), 30.25, tolerance=1e-7)
  expect_equal(env$extrinsic_drivers$evaluate("rainfall", -5.5), 30.25, tolerance=1e-7)
  
  ## interpolated range of points
  expect_equal(env$extrinsic_drivers$evaluate_range("rainfall", c(-7, 1, 7.8345)), c(49, 1, 61.37939025), tolerance=1e-6)
})
