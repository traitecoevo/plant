# Milestone C (#472 scope B / #537): CI-runnable validation of the scalar-
# templated FF16 demographic rate kernel (ff16_compute_rates_crown_top). The AD
# path is compiled into plant.so as forward-mode [[Rcpp::export]] free functions
# (header-only XAD, no reverse-mode tape, no on-the-fly Rcpp::sourceCpp), so this
# runs on CI like any other test. Two checks:
#   1. faithfulness -- the kernel reproduces the LIVE crown-centre fecundity rate
#      (so the AD result is a derivative of the real model, not a parallel formula);
#   2. gradient -- forward-mode d(fecundity_dt)/d(a_p1) matches a central finite
#      difference of the same kernel.
# The broader reverse-mode / emergent-output gradients (stand LAI, self-shading) are
# covered by test-ff16-stand-gradient.R and the gradient-regression fixture.

test_that("FF16 demographic rate kernel reproduces the live crown-centre rate", {
  ctrl <- Control(); ctrl$shading_model <- "crown-centre"
  s <- FF16_Strategy(); s$control <- ctrl
  ap1 <- s$pars$a_p1
  ind <- FF16_Individual(s)
  for (light_E in c(0.6, 0.9)) {
    env <- FF16_Environment(); env$set_fixed_environment(light_E, 1e4)
    for (height in c(8, 12, 15)) {
      ind$set_state("height", height)
      ind$compute_rates(env)
      kernel <- plant:::ff16_crown_top_fecundity_dt(height, light_E, ap1)
      # Bit-exact: the kernel is the single source the live path delegates to.
      expect_equal(kernel, ind$rate("fecundity"), tolerance = 1e-12)
    }
  }
})

test_that("forward-mode d(fecundity_dt)/d(a_p1) matches a finite difference", {
  ap1 <- FF16_Strategy()$pars$a_p1
  any_nonzero <- FALSE
  for (light_E in c(0.6, 0.9)) {
    for (height in c(8, 12, 15)) {
      ad <- plant:::ff16_fecundity_dt_grad_ap1(height, light_E)
      expect_true(is.finite(ad))
      e <- 1e-5 * ap1
      fd <- (plant:::ff16_crown_top_fecundity_dt(height, light_E, ap1 + e) -
             plant:::ff16_crown_top_fecundity_dt(height, light_E, ap1 - e)) / (2 * e)
      expect_equal(ad, fd, tolerance = 1e-6)
      if (abs(ad) > 0) any_nonzero <- TRUE
    }
  }
  # Guard against a vacuous pass (all gradients zero, e.g. if net production were
  # clamped everywhere): at least one configuration must exercise a real gradient.
  expect_true(any_nonzero)
})
