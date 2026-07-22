# R1 operator split of the TF24 soil water balance. The stiff per-layer drainage
# is integrated exactly (analytic_partial_flow); infiltration, the inter-layer
# cascade and uptake are the gentle residual (residual_rhs). Together they must
# reproduce the monolithic soil solve.

# monolithic soil rate, re-derived independently in R from the TF24 physics, as an
# uncontaminated reference: theta' = (water_in - drainage - uptake) / dz.
soil_defaults <- list(sat = 0.428, Ksat = 163.0411, npsi = 6.57, depth = 1.5,
                      a_infil = 1, b_infil = 8, theta_res = 1e-2, n = 5)

soil_K <- function(th, p) p$Ksat * (min(max(th, 0), p$sat) / p$sat)^(2 * p$npsi + 3)

mono_rate <- function(theta, rd, rain, p) {
  dz <- p$depth / p$n
  infil <- rain * max(0, 1 - p$a_infil * (theta[1] / p$sat)^p$b_infil)
  wout <- vapply(theta, soil_K, numeric(1), p = p)
  rate <- numeric(length(theta))
  for (i in seq_along(theta)) {
    win <- if (i == 1) infil else wout[i - 1]
    r <- (win - wout[i] - rd[i]) / dz
    if (theta[i] <= p$theta_res && !(r > 0)) r <- 0
    rate[i] <- r
  }
  rate
}

rk4 <- function(theta, rate_fn, dt, nsub) {
  h <- dt / nsub
  for (s in seq_len(nsub)) {
    k1 <- rate_fn(theta); k2 <- rate_fn(theta + 0.5 * h * k1)
    k3 <- rate_fn(theta + 0.5 * h * k2); k4 <- rate_fn(theta + h * k3)
    theta <- theta + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  theta
}

# The exact statement that the split IS the monolithic RHS: the drainage-flow rate
# (-K/dz) plus the residual equals the monolithic rate, everywhere in the interior.
testthat::test_that("the split decomposes the monolithic soil rate exactly", {
  p <- soil_defaults
  env <- Environment("TF24")
  env$set_soil_number_of_depths(p$n)
  dz <- p$depth / p$n
  cases <- list(
    list(theta = c(0.42, 0.38, 0.30, 0.22, 0.15), rd = c(3, 2, 1.5, 1, 0.5), rain = 0),
    list(theta = c(0.40, 0.35, 0.30, 0.25, 0.20), rd = rep(1, p$n),          rain = 20),
    list(theta = c(0.25, 0.24, 0.23, 0.22, 0.21), rd = c(2, 2, 2, 2, 2),     rain = 5)
  )
  for (cs in cases) {
    env$extrinsic_drivers_set_constant("rainfall", cs$rain)
    residual <- env$r_residual_rhs(cs$theta, cs$rd)
    flow_rate <- -vapply(cs$theta, soil_K, numeric(1), p = p) / dz
    mono <- mono_rate(cs$theta, cs$rd, cs$rain, p)
    expect_equal(residual + flow_rate, mono, tolerance = 1e-10)
  }
})

testthat::test_that("the exact drainage flow matches a fine integration", {
  p <- soil_defaults
  env <- Environment("TF24")
  env$set_soil_number_of_depths(p$n)
  dz <- p$depth / p$n
  for (th0 in c(0.42, 0.30, 0.15)) {
    for (dt in c(0.5, 2.0)) {
      exact <- env$r_analytic_partial_flow(rep(th0, p$n), dt)
      fine <- rk4(rep(th0, p$n), function(th) -vapply(th, soil_K, numeric(1), p = p) / dz,
                  dt, nsub = 20000)
      expect_lt(max(abs(exact - fine)), 1e-7)
      expect_true(all(exact > 0))                       # positivity by construction
    }
  }
})

testthat::test_that("the touchdown time brings a layer to the residual floor", {
  p <- soil_defaults
  env <- Environment("TF24")
  env$set_soil_number_of_depths(p$n)
  for (th0 in c(0.42, 0.25)) {
    t_touch <- env$r_drainage_touchdown_time(th0, 0)
    expect_true(is.finite(t_touch) && t_touch > 0)
    at_floor <- env$r_analytic_partial_flow(rep(th0, p$n), t_touch)
    expect_equal(at_floor[1], p$theta_res, tolerance = 1e-6)
  }
  expect_true(is.infinite(env$r_drainage_touchdown_time(soil_defaults$theta_res, 0)))
})

# Strang integration converges to the monolithic solve. Where the monolithic
# explicit solve is well resolved (drought, mild forcing, at the floor) the split
# matches it closely; in the monsoon regime the explicit reference is itself
# stability-limited (the reason R1 exists), so there we only require the split to
# converge as the step shrinks and to stay bounded and positive.
split_run <- function(env, theta0, rd, T_end, nstep) {
  theta <- theta0
  h <- T_end / nstep
  for (step in seq_len(nstep)) {
    theta <- env$r_analytic_partial_flow(theta, h / 2)
    theta <- rk4(theta, function(th) env$r_residual_rhs(th, rd), h, nsub = 4)
    theta <- env$r_analytic_partial_flow(theta, h / 2)
  }
  theta
}

testthat::test_that("Strang split matches the monolithic solve when it is resolved", {
  p <- soil_defaults
  env <- Environment("TF24")
  env$set_soil_number_of_depths(p$n)
  T_end <- 4
  scenarios <- list(
    list(rain = 0, rd = rep(2.0, p$n),        theta0 = rep(0.40, p$n)),  # drought
    list(rain = 1, rd = c(3, 2.5, 2, 1.5, 1), theta0 = rep(0.30, p$n)),  # mild
    list(rain = 5, rd = rep(6.0, p$n),        theta0 = rep(0.08, p$n))   # near the floor
  )
  for (s in scenarios) {
    env$extrinsic_drivers_set_constant("rainfall", s$rain)
    mono <- rk4(s$theta0, function(th) mono_rate(th, s$rd, s$rain, p), T_end, nsub = 16000)
    e_coarse <- max(abs(split_run(env, s$theta0, s$rd, T_end, 100) - mono))
    e_fine   <- max(abs(split_run(env, s$theta0, s$rd, T_end, 400) - mono))
    expect_lt(e_fine, e_coarse)          # Strang splitting error -> 0 with h
    expect_lt(e_fine, 1e-2)              # and is small at a modest step
  }
})

testthat::test_that("Strang split converges and stays bounded in the stiff monsoon", {
  p <- soil_defaults
  env <- Environment("TF24")
  env$set_soil_number_of_depths(p$n)
  env$extrinsic_drivers_set_constant("rainfall", 30)
  T_end <- 4; rd <- rep(1.0, p$n); theta0 <- rep(0.20, p$n)
  mono <- rk4(theta0, function(th) mono_rate(th, rd, 30, p), T_end, nsub = 16000)
  e1 <- max(abs(split_run(env, theta0, rd, T_end, 400) - mono))
  e2 <- max(abs(split_run(env, theta0, rd, T_end, 1600) - mono))
  expect_lt(e2, e1)                                    # converging toward the reference
  fine <- split_run(env, theta0, rd, T_end, 1600)
  expect_true(all(is.finite(fine) & fine > 0 & fine < p$sat + 0.05))
})
