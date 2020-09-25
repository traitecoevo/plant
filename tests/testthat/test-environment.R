
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  e <- environment_types[[x]]

  context(sprintf("Environment-%s",x))

  test_that("Empty environment", {
    p <- Parameters(x, e)()
    env <- make_environment(x, p)

    ## At this point, we should have full canopy openness, partly because
    ## the spline is just not constructed.
    expect_equal(env$canopy_openness(0), 1.0)
    expect_equal(env$canopy_openness(100), 1.0)

    spline <- env$canopy$canopy_interpolator
    expect_equal(spline$size, 33)
    expect_equal(spline$x, seq(0,1, length.out=33))
  })

  test_that("Manually set environment", {
    env <- make_environment(x, Parameters(x, e)())
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

  test_that("Disturbance related parameters", {
    env <- make_environment(x, Parameters(x, e)())
    expect_identical(env$time, 0.0)
    expect_identical(env$patch_survival_conditional(env$time), 1.0)

    disturbance <- Disturbance(30.0)
    env$time <- 10
    expect_identical(env$patch_survival_conditional(0), disturbance$pr_survival_conditional(env$time, 0))
    expect_identical(env$patch_survival_conditional(2), disturbance$pr_survival_conditional(env$time, 2))

    expect_is(env$disturbance_regime, "Disturbance")
  })

  test_that("Seed rain related parameters", {
    env <- make_environment(x, Parameters(x, e)())
    expect_error(env$seed_rain_dt, "Cannot get seed rain for empty environment")

    z <- c(.1, .2)
    env <- test_environment(x, 10, n_strategies=2, seed_rain=z)

    expect_identical(env$seed_rain_dt, z[[1]])
    env$set_seed_rain_index(1)
    expect_identical(env$seed_rain_dt, z[[1]])
    env$set_seed_rain_index(2)
    expect_identical(env$seed_rain_dt, z[[2]])
    expect_error(env$set_seed_rain_index(0), "Invalid value for index")
    expect_error(env$set_seed_rain_index(3), "Index 3 out of bounds")
  })
}
