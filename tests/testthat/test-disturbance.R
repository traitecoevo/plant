context("Disturbance")

test_that("Base class", {
  obj <- Disturbance_Regime()
  expect_is(obj, "Disturbance_Regime")
  expect_is(obj, "R6")

  expect_error(obj$density(), "argument \"time\" is missing, with no default")
  expect_true(is.na(obj$density(1)))
  expect_true(is.na(obj$pr_survival(1)))
})

test_that("No disturbance", {
  obj <- NoDisturbance()
  
  # always 1 for individual patches
  expect_equal(obj$pr_survival(1), 1.0)
  expect_equal(obj$pr_survival(100), 1.0)

  # integration test to show constant density  
  time <- seq(0, 100, len = 1e3)
  d <- obj$density(time)
  expect_equal(trapezium(time, d), 100)
})

test_that("Weibull disturbance regime", {
  # max_patch_lifetime determines the shape of the disturbance distribution
  m <- 105.32
  obj <- WeibullDisturbance_Regime(m)
  
  # check mean interval
  expect_equal(obj$mean_interval(), 30)
  
  # A Weibull distribution has PDF proportional to
  #   (x/lambda)^(k-1) * exp(-(x/lambda)^k)
  # however, we use the reparameterisation
  #   b = lambda^(-k)
  #   bkx^(k-1) * exp(-bx^k)
  k = 2
  lambda = 33.85138
  b = lambda^(-k)
  
  # that said, we use a scaled inverse CDF for our patch density
  expect_false(obj$density(m) == dweibull(m, k, lambda))
  
  p0 <- k * (b^(1.0 / k) / gamma(1.0 / k))
  expect_equal(obj$density(m), p0 * (1 - pweibull(m, k, lambda)))
  
  # integration test to show normalised density  
  time <- seq(0, 100, len = 1e3)
  d <- obj$density(time)
  expect_equal(trapezium(time, d), 0.9999, tolerance = 0.0001)
  
  # check cdf
  expect_equal(obj$cdf(m), 0.9999375, tolerance = 0.0001)  
  expect_equal(obj$cdf(m), pweibull(m, k, lambda))
  
  # check icdf - internal calibration chosen to match Daniels prev. implementation
  expect_equal(obj$icdf(1 - 0.9999375), m, tolerance = 0.0001)  
  expect_equal(obj$icdf(6.25302620663814e-05), m)  
  })

# TODO: check what this is for
# test_that("Reference survival eps gives correct running time", {
#   ## From falster-traitdiversity: src/base/ebt/site.cpp,
#   ## site::solve_patchage_dist():
#   f1 <- function(mean) {
#     2.633 * mean / 3.0 * 4.0
#   }
#
#   ## Wrapper for getting same out of our disturbance class:
#   ctrl <- Control()
#   eps <- ctrl$schedule_patch_survival
#   f2 <- function(mean) {
#     Disturbance(mean)$cdf(eps)
#   }
#
#   age <- seq(1, 200, length.out=101)
#   y1 <- f1(age)
#   y2 <- sapply(age, f2)
#
#   expect_equal(y1, y2)
# })
