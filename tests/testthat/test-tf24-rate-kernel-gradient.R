# Phase F1-full (#472 scope B): CI-runnable validation of the scalar-templated
# TF24 demographic rate kernel (tf24_compute_rates_from_net). The AD path is
# compiled into plant.so as forward-mode [[Rcpp::export]] free functions
# (header-only XAD, no reverse-mode tape, no on-the-fly Rcpp::sourceCpp), so this
# runs on CI like any other test. Two checks:
#   1. faithfulness -- the kernel rate fill reproduces the LIVE crown-centre
#      fecundity rate (so the AD result is a derivative of the real model, not a
#      parallel formula); TF24_Strategy's rate methods delegate to the kernel;
#   2. gradient -- forward-mode d(fecundity_dt)/d(vcmax_25), with the leaf-profit
#      sensitivity (Leaf::dprofit_dvcmax25) injected into net, matches a central
#      finite difference of the live crown-centre net through the same kernel.
# The broader reverse-mode 27-trait sweep + emergent SCM gradient are covered by
# test-tf24-offspring-gradient.R and test-tf24f-census-gradient.R.

test_that("TF24 demographic rate kernel reproduces the live crown-centre rate", {
  ctrl <- Control(); ctrl$shading_model <- "crown-centre"
  ctrl$GSS_tol_abs <- 1e-9   # match the kernel free function's tight collar optimum
  s <- TF24_Strategy(); s$control <- ctrl
  vcmax <- s$pars$vcmax_25
  ind <- TF24_Individual(s)
  for (light_E in c(0.6, 0.9)) {
    env <- TF24_Environment(); env$set_fixed_environment(light_E, 1e4)
    for (height in c(8, 12, 15)) {
      ind$set_state("height", height)
      ind$compute_rates(env)
      kernel <- plant:::tf24_crown_centre_fecundity_dt(height, light_E, vcmax)
      # Bit-exact: the kernel is the single source the live path delegates to,
      # given the same (live, crown-centre) net production.
      expect_equal(kernel, ind$rate("fecundity"), tolerance = 1e-12)
    }
  }
})

test_that("forward-mode d(fecundity_dt)/d(vcmax_25) matches a finite difference", {
  vcmax <- TF24_Strategy()$pars$vcmax_25
  any_nonzero <- FALSE
  # Heights where reproduction is active (near/above hmat=16.6): below that the
  # fecundity rate -- and its gradient -- is ~0, so a relative check is vacuous.
  # vcmax_25 flows through the hydraulic leaf optimisation, so the FD floor is the
  # leaf root-find noise (~1e-8 at the matched step), not the closed-form ~1e-12.
  for (light_E in c(0.6, 0.9)) {
    for (height in c(12, 15)) {
      ad <- plant:::tf24_fecundity_dt_grad_vcmax(height, light_E)
      expect_true(is.finite(ad))
      e <- 1e-6 * vcmax
      fd <- (plant:::tf24_crown_centre_fecundity_dt(height, light_E, vcmax + e) -
             plant:::tf24_crown_centre_fecundity_dt(height, light_E, vcmax - e)) / (2 * e)
      # Portable tolerance: the FD reference goes through the hydraulic leaf root-find,
      # whose noise floor is compiler/libm-dependent (~7e-5 on some platforms vs ~1e-8 on
      # the dev machine), so an AD-vs-FD check must be loose enough to be cross-platform.
      expect_equal(ad, fd, tolerance = 1e-3)
      if (abs(ad) > 0) any_nonzero <- TRUE
    }
  }
  # Guard against a vacuous pass (all gradients zero, e.g. if net production were
  # clamped everywhere): at least one configuration must exercise a real gradient.
  expect_true(any_nonzero)
})
