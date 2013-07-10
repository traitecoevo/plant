source("helper-tree.R")

context("Environment")

e <- new(Environment, new(Parameters))

## At this point, we should have full canopy openness, partly because
## the spline is just not constructed.
expect_that(e$canopy_openness(0), equals(1.0))
expect_that(e$canopy_openness(100), equals(1.0))

spline <- e$light_environment
expect_that(spline$size, equals(0))
expect_that(spline$x, equals(numeric(0)))

## Now, set the light environment.
hh <- seq(0, 10, length=101)
light.env <- function(x)
  exp(x/(max(hh)*2)) - 1 + (1 - (exp(.5) - 1))/2
ee <- light.env(hh)
env <- new(Spline)
env$init(hh, ee)

## And set it
e$light_environment <- env

expect_that(e$light_environment$xy, is_identical_to(env$xy))

hmid <- (hh[-1] + hh[-length(hh)])/2
expect_that(sapply(hmid, e$light_environment$eval),
            is_identical_to(sapply(hmid, env$eval)))

expect_that(e$time, is_identical_to(0.0))
expect_that(e$patch_survival(e$time), is_identical_to(1.0))

disturbance <- new(Disturbance, 30.0)
e$time <- 10
expect_that(e$patch_survival(0),
            is_identical_to(disturbance$survival_probability(0, e$time)))
expect_that(e$patch_survival(2),
            is_identical_to(disturbance$survival_probability(2, e$time)))


rain.blank <- e$seed_rain
expect_that(rain.blank$get(), throws_error())
expect_that(e$seed_rain_rate, throws_error())

x <- c(.1, .2)
e <- test.environment(10, n.strategies=2, seed.rain=x)
expect_that(e$seed_rain, is_identical_to(x))
expect_that(e$seed_rain_rate, is_identical_to(x[[1]]))
e$set_seed_rain_index(1)
expect_that(e$seed_rain_rate, is_identical_to(x[[1]]))
e$set_seed_rain_index(2)
expect_that(e$seed_rain_rate, is_identical_to(x[[2]]))
expect_that(e$set_seed_rain_index(0), throws_error())
expect_that(e$set_seed_rain_index(3), throws_error())
