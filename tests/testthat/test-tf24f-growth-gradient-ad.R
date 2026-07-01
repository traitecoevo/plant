# Exact-AD growth-rate gradient for TF24f (#472 scope B / #537 A1). TF24f overrides
# Strategy::growth_rate_gradient_height_ad with an analytic d(dheight/dt)/d(height) at
# the FIXED tracked collar (forward-AD of the demographic growth kernel, with profit
# carrying its height derivative through the kmax + light + E_up channels). This is the
# direct replacement for Node::growth_rate_gradient's 1e-6 backward FD, reachable in
# run_scm via control(node_gradient_exact_ad = TRUE) -- now exposed to R (the C++ field
# existed since #537 A1 but was never in the Control interface, so neither this nor
# FF16's exact gradient could be switched on from run_scm before).
#
# SCOPE. Validated in a fixed env here (matches a fine FD of the live growth rate to
# ~1e-6). It is also correct on the real crown-centre RESIDENT light: an in-run_scm
# (native-env) check found the analytic dprofit_dh matches a native profit-FD to ~1e-4
# and the AD g' matches a native height_dt-FD to ~1e-5 even for the tallest cohort. (An
# earlier "~18% off" reading was an artifact of validating through an Rcpp::as<>
# round-tripped environment, which is NOT faithful for the crown-sampled light above
# cohort heights -- not a fault in the gradient.) node_gradient_exact_ad stays OFF by
# default for bit-compatibility; the exact-AD
# gradient does NOT change the TF24f census fidelity floor (that is a separate,
# Rcpp::as<>-limited reconstruction issue at long horizons, not g').

test_that("TF24f growth_rate_gradient_exact matches a fine FD in a fixed env", {
  s <- TF24f_Strategy()
  ind <- TF24f_Individual(s)
  env <- Environment("TF24f")
  env$set_fixed_environment(0.85, 1e4)
  collar <- 0.9   # a representative tracked root-collar potential (MPa magnitude)

  growth_rate <- function(h) {
    ind$set_state("height", h)
    ind$set_state("opt_root_psi_state", collar)
    ind$compute_rates(env)
    ind$rate("height")
  }

  for (height in c(2, 5, 9)) {
    ind$set_state("height", height)
    ind$set_state("opt_root_psi_state", collar)
    ind$compute_rates(env)
    g <- ind$growth_rate_gradient_exact(env)
    expect_true(is.finite(g))
    e <- 1e-3 * height
    fd <- (growth_rate(height + e) - growth_rate(height - e)) / (2 * e)
    expect_equal(g, fd, tolerance = 1e-5)
  }
})

test_that("growth_rate_gradient_exact is NA for TF24 (no AD override)", {
  # TF24 (and K93) inherit the Strategy<E> base default (NA), so Node falls back to
  # the finite difference; only TF24f and FF16 provide the exact gradient.
  ind <- TF24_Individual(TF24_Strategy())
  env <- Environment("TF24")
  env$set_fixed_environment(0.85, 1e4)
  ind$set_state("height", 5)
  ind$compute_rates(env)
  expect_true(is.na(ind$growth_rate_gradient_exact(env)))
})

test_that("node_gradient_exact_ad is exposed and run_scm honours it for TF24f", {
  # The Control field is now in the R interface; a TF24f patch runs end-to-end with
  # the exact-AD growth-rate gradient enabled (kept short so it is CI-cheap).
  mk <- function(exact_ad) {
    p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- 4L
    p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                        birth_rate = list(20))
    p$node_schedule_times <- list(seq(0, 4L, length.out = 9L))
    ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                    ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE,
                    node_gradient_exact_ad = exact_ad)
    run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)
  }
  expect_identical(Control(node_gradient_exact_ad = TRUE)$node_gradient_exact_ad, TRUE)
  fd <- mk(FALSE); ad <- mk(TRUE)
  # Both produce a finite, positive stand; at this short horizon (smooth interior
  # cohorts) the exact-AD log-density tracks the FD closely.
  expect_true(is.finite(ad$offspring_production[[1]]) && ad$offspring_production[[1]] > 0)
  ldf <- fd$patch$species[[1]]$log_densities
  lda <- ad$patch$species[[1]]$log_densities
  expect_lt(max(abs(lda - ldf)), 1e-2)
})

test_that("TF24f AD acclimation gradient is drought-robust on long patches (#527)", {
  # Regression for the lifetime>=8 abort: TF24f's AD acclimation gradient
  # (solve_leaf, default use_ad_gradient = TRUE) used to call
  # Leaf::dprofit_droot_collar_psi at the drought-SHUTDOWN operating point
  # (collar = -psi_crit), where the non-extrapolating transport spline's
  # derivative is out of domain -- a hard "Extrapolation disabled" abort on long /
  # dry patches. solve_leaf now gates on prepare_collar_solve (as the #526 FD path
  # does) and returns a 0 gradient at shutdown, so a long TF24f patch runs to
  # completion with the exact gradient, agreeing with the FD-gradient run.
  mk <- function(H, use_ad) {
    p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
    p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                        birth_rate = list(20))
    s <- p$strategies[[1]]; s$use_ad_gradient <- use_ad; p$strategies[[1]] <- s
    p$node_schedule_times <- list(seq(0, H, length.out = 9L))
    ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                    ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
    run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)
  }
  for (H in c(8L, 12L)) {
    ad <- mk(H, TRUE)               # the path that used to abort
    expect_true(is.finite(ad$offspring_production[[1]]) &&
                ad$offspring_production[[1]] > 0)
    fd <- mk(H, FALSE)              # the FD path that never aborted
    expect_equal(ad$offspring_production[[1]], fd$offspring_production[[1]],
                 tolerance = 5e-2)
  }
})
