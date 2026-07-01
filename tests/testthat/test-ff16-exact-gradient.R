# Milestone C (#472 scope B / #537 A1): the exact AD growth-rate gradient as a
# first-class R method. Individual$growth_rate_gradient_exact(env) returns the
# strategy's exact d(growth rate)/d(height) via forward-mode AD (no on-the-fly
# compilation -- it is compiled into the package), so this runs on CI. Validated
# against a fine finite difference of the real growth rate (set height ->
# compute_rates -> rate("height")). This is the quantity Node::growth_rate_gradient
# obtains by FD, available exactly when control$node_gradient_exact_ad is set.

testthat::test_that("Individual$growth_rate_gradient_exact matches a fine FD (FF16, A1)", {
  s <- FF16_Strategy()
  ind <- FF16_Individual(s)
  env <- FF16_Environment()
  env$set_fixed_environment(0.85, 1e4)

  growth_rate <- function(h) {
    ind$set_state("height", h)
    ind$compute_rates(env)
    ind$rate("height")
  }

  for (height in c(4, 8, 12)) {
    ind$set_state("height", height)
    g <- ind$growth_rate_gradient_exact(env)
    expect_true(is.finite(g))
    e <- 1e-6
    fd <- (growth_rate(height + e) - growth_rate(height - e)) / (2 * e)
    expect_equal(g, fd, tolerance = 1e-5)
  }
})

testthat::test_that("growth_rate_gradient_exact is NA for strategies without an AD gradient", {
  # K93 / TF24 inherit the Strategy<E> base default (NA), so Node falls back to FD.
  ind <- K93_Individual(K93_Strategy())
  env <- K93_Environment()
  env$set_fixed_environment(0.5, 1e4)
  ind$set_state("height", 5)
  expect_true(is.na(ind$growth_rate_gradient_exact(env)))
})
