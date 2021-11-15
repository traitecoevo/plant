
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
  context(sprintf("Rainfall-FF16",x))
  env <- make_environment("FF16")
  
  ## simple quadratic
  x <- seq(-10, 10, 2)
  y <- x^2
  env$rainfall_init(x, y)
  
  ## discrete points used for spline creation
  expect_equal(env$rainfall_eval(2), 4)
  expect_equal(env$rainfall_eval(-2), 4)
  
  ## interpolated points
  expect_equal(env$rainfall_eval(3), 9, tolerance=1e-2)
  expect_equal(env$rainfall_eval(-3), 9, tolerance=1e-2)
  expect_equal(env$rainfall_eval(5.5), 30.25, tolerance=1e-2)
  expect_equal(env$rainfall_eval(-5.5), 30.25, tolerance=1e-2)
  
  ## interpolated range of points
  expect_equal(env$rainfall_eval_range(c(-7, 1, 7.8345)), c(49, 1, 61.37939025), tolerance=1e-2)
  
  ## clear points, set new points and recompute spline
  env$rainfall_clear_points()
  x <- seq(-10, 10, 1)
  y <- 2*x
  env$rainfall_add_points(x, y)
  env$rainfall_recompute()
  expect_equal(env$rainfall_eval(3), 6, tolerance=1e-2)
  expect_equal(env$rainfall_eval(3.3), 6.6, tolerance=1e-2)
})
