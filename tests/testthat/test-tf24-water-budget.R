# Water-budget closure for the TF24 multi-layer soil column, and the
# non-negativity floor on the rainfall driver.
#
# Coverage before this file existed was a single-layer closure check inside
# test-environment-TF24.R that omitted root uptake (it ran for 0.01 yr, so
# uptake was negligible). These tests close the budget with uptake included,
# over several layer counts and across the wet/dry range -- including at and
# below the residual-moisture floor, where the guard in
# TF24_Environment::compute_rates could plausibly have leaked mass.

# Total column water (m) per output step, and the closure residual.
#
# The budget is a storage form that avoids having to reconstruct runoff:
#   storage_start + infiltration == storage_end + drainage + uptake
# Runoff never enters the column, so it cancels; `sum_rainfall` is only needed
# to check the driver itself (see the rainfall-floor tests below).
tf24_water_budget <- function(n_layers = 5, theta_0 = 0.25, rainfall = 1.0,
                              lifetime = 2, lma = 0.0825) {
  env <- Environment("TF24")
  env$set_soil_number_of_depths(n_layers)
  env$set_soil_water_state(rep(theta_0, n_layers))
  env$extrinsic_drivers_set_constant("rainfall", rainfall)

  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(lma, "lma"))

  out <- run_scm(p, env, collect = TRUE)

  # Layers are uniform, so the first cumulative depth is the layer thickness.
  dz <- out$env$soil_depth$soil_depth[[1]]
  moist <- out$env$soil_moist
  storage <- tapply(moist$soil_moist * dz, moist$step, sum)

  flux <- out$env$soil_moist_cumulative_flux
  final <- flux[nrow(flux), ]

  start <- storage[[as.character(min(flux$step))]]
  end   <- storage[[as.character(max(flux$step))]]

  supplied <- start + final$sum_infiltration
  accounted <- end + final$sum_drainage + final$sum_resource_depletion

  list(
    supplied = supplied,
    accounted = accounted,
    residual = supplied - accounted,
    rel_residual = (supplied - accounted) / supplied,
    min_theta = min(moist$soil_moist),
    max_theta = max(moist$soil_moist),
    n_layers = n_layers
  )
}

test_that("TF24 water budget closes across layer counts", {
  # 1 layer is the previously-covered case; 5 is the constructor default; 15
  # matches the layer count LPJ-GUESS and comparable land models use over 1.5 m.
  #
  # The per-layer `rel_residual` check is the discretisation-independence claim
  # too: a change in discretisation must not create or destroy water, so the
  # residual has to stay at round-off rather than grow with the layer count.
  # Asserting it inside this loop covers that with one set of runs; a separate
  # test re-running the same three configurations to re-assert the same bound
  # cost ~2.8 s for no additional coverage.
  for (n in c(1, 5, 15)) {
    b <- tf24_water_budget(n_layers = n)
    expect_equal(b$supplied, b$accounted,
                 info = sprintf("n_layers = %d", n))
    expect_lt(abs(b$rel_residual), 1e-12)
  }
})

test_that("TF24 water budget closes over the wet-to-dry range", {
  cases <- list(
    list(theta_0 = 0.428, rainfall = 2.0, label = "saturated, heavy rain"),
    list(theta_0 = 0.250, rainfall = 1.0, label = "moist, default rain"),
    list(theta_0 = 0.050, rainfall = 0.0, label = "dry, no rain")
  )
  for (case in cases) {
    b <- tf24_water_budget(theta_0 = case$theta_0, rainfall = case$rainfall,
                           lifetime = 5)
    expect_equal(b$supplied, b$accounted, info = case$label)
  }
})

test_that("the residual-moisture guard does not leak water", {
  # TF24_Environment::compute_rates refuses to dry a layer already at or below
  # soil_moist_residual (1e-2). That clamp discards a rate without accounting
  # for it, so it *could* break closure. These start the column at and below the
  # floor to exercise it directly.
  for (theta_0 in c(0.010, 0.005)) {
    b <- tf24_water_budget(theta_0 = theta_0, rainfall = 0.0, lifetime = 10)
    expect_equal(b$supplied, b$accounted,
                 info = sprintf("theta_0 = %g", theta_0))
  }
})

test_that("Campbell transport shuts off before the residual floor is reached", {
  # Documents *why* the guard above never leaks in practice: with n_psi = 6.57
  # the conductivity exponent is 2*n_psi+3 = 16.14, so K(theta) collapses and
  # psi(theta) diverges far above theta_r. A column started at theta = 0.05 with
  # no rain neither drains nor is drawn down measurably -- transport has already
  # stopped, so the clamp is a safety net rather than an active mass sink.
  env <- TF24_Environment()
  # K(theta) = K_sat * (theta/theta_sat)^(2*n_psi+3); soil_K_from_soil_theta is
  # not exposed to R, so evaluate the same expression from the exposed pars.
  k_at <- function(theta) {
    env$K_sat * (theta / env$soil_moist_sat)^(2 * env$n_psi + 3)
  }
  expect_lt(k_at(0.05), 1e-10)                    # m yr^-1, against K_sat = 163
  # psi at theta = 0.05 is far past any plausible root shutoff (~6 MPa).
  expect_gt(env$psi_from_soil_moist(0.05), 100)   # MPa

  b <- tf24_water_budget(theta_0 = 0.05, rainfall = 0.0, lifetime = 20)
  expect_equal(b$min_theta, 0.05, tolerance = 1e-6)
})

test_that("rainfall driver is floored at zero", {
  # Extrinsic drivers are interpolated with a cubic spline, which undershoots
  # badly on intermittent forcing. Unfloored, a negative rainfall value gives
  # negative infiltration and an unphysical drying rate; where the residual
  # guard then clamps that rate, the removal is recorded in sum_rainfall but
  # never applied and the budget stops closing.
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.25, 5))
  # A driver that is explicitly negative over part of its range.
  times <- seq(0, 2, length.out = 21)
  env$extrinsic_drivers_set_variable("rainfall", times, rep(c(-1, 1), length.out = 21))
  expect_lt(min(env$extrinsic_drivers_evaluate_range("rainfall", seq(0, 2, 0.01))), 0)

  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- 2
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  out <- run_scm(p, env, collect = TRUE)

  flux <- out$env$soil_moist_cumulative_flux
  # Cumulative rainfall and infiltration are both monotone non-decreasing: the
  # column is never charged a negative rainfall or a negative infiltration.
  expect_true(all(diff(flux$sum_rainfall) >= 0))
  expect_true(all(diff(flux$sum_infiltration) >= 0))

  # And the budget still closes under a driver that dips negative.
  dz <- out$env$soil_depth$soil_depth[[1]]
  moist <- out$env$soil_moist
  storage <- tapply(moist$soil_moist * dz, moist$step, sum)
  final <- flux[nrow(flux), ]
  expect_equal(storage[[as.character(min(flux$step))]] + final$sum_infiltration,
               storage[[as.character(max(flux$step))]] + final$sum_drainage +
                 final$sum_resource_depletion)
})

test_that("an interpolated driver cannot undershoot its own series", {
  env <- TF24_Environment()
  # Intermittent daily rainfall: the case that motivated the floor, and the case
  # this assertion used to record as broken. A driver arrives as values with no
  # slopes, and the rule that chooses them decides whether a wet day can push the
  # days beside it below zero. It can no longer: the fit is monotone on every span,
  # so a read is bounded by the two control points that bracket it.
  #
  # What that replaced, on this exact series: 12% of points negative, reaching
  # -1.86 mm, with the integral conserved -- so a total-rainfall check could not
  # see it and this diagnostic had to exist to surface it. The diagnostic stays,
  # because a caller can still supply a series that is itself negative.
  set.seed(42)
  n <- 365 * 3
  times <- seq(0, 3, length.out = n)
  doy <- times %% 1
  y <- ifelse(doy > 0.5 & doy < 0.75 & stats::runif(n) < 0.35,
              stats::rexp(n, 1 / 8), 0)
  env$extrinsic_drivers_set_variable("rainfall", times, y)

  expect_silent(
    res <- suppressMessages(check_driver_interpolation(env, "rainfall", times, y)))

  expect_equal(min(y), 0)
  expect_gte(res$min, 0)
  expect_equal(res$frac_negative, 0)
  expect_equal(res$negative_area, 0)
  # And it stays a faithful total, which the overshooting fit also managed -- by
  # cancelling one error against another rather than by making none.
  expect_equal(res$integral_evaluated, res$integral_supplied, tolerance = 1e-5)
  expect_lte(res$max, max(y))
})

test_that("check_driver_interpolation is quiet on a smooth series", {
  env <- TF24_Environment()
  times <- seq(0, 5, length.out = 120)
  y <- 0.4 + 0.4 * sin(2 * pi * times)
  y[y < 0] <- 0
  env$extrinsic_drivers_set_variable("rainfall", times, y)
  expect_silent(suppressMessages(check_driver_interpolation(env, "rainfall", times, y)))
})
