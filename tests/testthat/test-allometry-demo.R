# Guards the numbers TF24_flexible_allometry_demo.qmd asserts in prose.
#
# The demo's chunks all run, but a chunk that runs is not a chunk that is
# checked: prose beside a number survives the number being falsified. One claim
# in the first draft did exactly that (a Huber-value rise quoted from a shorter
# drought than the one plotted), which is why this file exists.

demo_helpers <- test_path("..", "..", "overstorey_staging",
                          "allometry_demo_helpers.R")

test_that("the flexible-allometry demo's claims still hold", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  history <- demo_soil_history(wet = 0.30, dry = 0.11,
                               years_wet = 2, years_dry = 4)
  fixed <- allometry_trajectory(0.0, history, dt = 0.25)
  flex  <- allometry_trajectory(1.0, history, dt = 0.25)

  ## Tolerances are loose because these are integrated trajectories, not
  ## closed forms; what is being guarded is the claim, not the digits.

  ## "thins to about 19 per cent of what its height prefers"
  expect_equal(min(flex$canopy_fraction), 0.19, tolerance = 0.15)
  ## the fixed plant does not move off its allometry at all
  expect_equal(max(abs(fixed$canopy_fraction - 1)), 0, tolerance = 1e-8)

  ## "rises to about 2.4 times its preferred value"
  expect_equal(max(flex$huber_ratio), 2.4, tolerance = 0.15)
  ## and psi >= 0 is invariant, so the ratio never drops below 1
  expect_gte(min(flex$huber_ratio), 1 - 1e-8)

  ## "rebuilds to ~98 per cent of preferred"
  expect_gt(tail(flex$canopy_fraction, 1), 0.9)

  ## Height never falls, in either model.
  expect_true(all(diff(flex$height) >= -1e-10))
  expect_true(all(diff(fixed$height) >= -1e-10))

  ## THE ISSUE THE DEMO EXISTS TO FLAG: thinning buys ~1 per cent survival and
  ## costs ~2 m of height. If either of these moves materially, the demo's
  ## conclusion has changed and its prose needs rewriting -- which is the point
  ## of asserting them.
  expect_equal(max(flex$survivorship / fixed$survivorship), 1.01,
               tolerance = 0.05)
  expect_equal(tail(fixed$height, 1) - tail(flex$height, 1), 2.07,
               tolerance = 0.4)
})

test_that("the demo's gate map brackets the replacement gate", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  gate <- allometry_gate_map(a_pl0 = 1.0, height = 5, n = 21L)

  ## Opposite signs: leaf area thins while sapwood per leaf area rises.
  expect_lt(gate$dphi[[1]], 0)
  expect_gt(gate$dpsi[[1]], 0)

  ## "by r = 0.33 the rates are already ~1e-7" -- the gate is sharp, which is
  ## the demo's explanation for why the response is small in most stands.
  near_third <- gate[which.min(abs(gate$r - 1 / 3)), ]
  expect_lt(abs(near_third$dphi), 1e-5)

  ## Ample reserves: no thinning at all.
  expect_equal(gate$dphi[[nrow(gate)]], 0, tolerance = 1e-12)

  ## Height growth is non-negative throughout.
  expect_true(all(gate$dheight >= 0))
})

# An independent R implementation of the equations written out in the demo's
# "The mathematics" section, checked against what the C++ actually computes.
#
# This is the guard that matters most for a documented mechanism: the demo can
# only be trusted if its equations ARE the model. Written from the prose rather
# than from the source, so a divergence between the two shows up here.
#
# Constants not exposed to R are repeated with the member they mirror named. If
# one of those changes, this test fails -- which is the intent, not a nuisance:
# the prose would need changing too.
tf24_departure_rates_r <- function(s, height, phi, psi, storage, P, A) {
  p <- s$pars
  storage_gate_width      <- 0.1    # TF24_Strategy::storage_gate_width
  storage_prod_eps        <- 1e-4   # TF24_Strategy::storage_prod_eps
  plasticity_gate_width   <- 0.02   # TF24_Strategy::plasticity_gate_width

  eta_c <- 1 - 2 / (1 + p$eta) + 1 / (1 + 2 * p$eta)   # CanopyShape::eta_c

  ## Sizes
  A_s   <- p$theta * A * exp(psi)
  m_s   <- A_s * height * eta_c * p$rho
  S_max <- p$a_st1 * m_s
  r     <- if (S_max > 0) storage / S_max else 0

  ## The replacement gate, and the sapwood share scaled by how over-built it is
  rho_l <- 1 - p$a_pl0 / (1 + exp((r - p$a_pl1) / plasticity_gate_width))
  rho_s <- rho_l * exp(-psi / p$a_pl2)
  u_l   <- 1 - rho_l
  u_s   <- 1 - rho_s

  ## Growth flux, and its split
  G     <- 1 / (1 + exp(-(r - p$a_st2) / storage_gate_width))
  Ppos  <- 0.5 * (P + sqrt(P^2 + storage_prod_eps^2))
  F     <- Ppos * G
  f_r   <- p$a_f1 / (1 + exp(p$a_f2 * (1 - height / p$hmat)))
  f_g   <- 1 - f_r
  sigma <- rho_l * (1 - exp(phi / p$a_pl2))

  ## Rebuilding buys the whole package at the preferred ratio
  c_reb <- p$lma + p$a_r1 + p$theta * height * eta_c * p$rho * (1 + p$a_b1)
  b     <- sigma * F * f_g / (c_reb * A)

  list(dphi = b - u_l * p$k_l,
       dpsi = u_l * p$k_l - u_s * p$k_s + b * (exp(-psi) - 1),
       r = r, rho_l = rho_l, rho_s = rho_s, b = b)
}

test_that("the demo's equations reproduce the C++ departure rates", {
  skip_on_cran()
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.2, 5))
  env$set_fixed_environment(1.0, 40)

  s <- TF24_Strategy(collect_all_auxiliary = TRUE)
  s$pars$a_pl0 <- 1.0

  ## A grid that exercises both departures, both signs of production, and the
  ## whole range of the gate -- including the resting point and states far off
  ## the trajectory.
  grid <- expand.grid(height = c(1, 5, 15),
                      phi = c(0, -0.4, -1.5),
                      psi = c(0, 0.3, 1.2),
                      storage_frac = c(0, 0.03, 0.5))

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    ind <- TF24_Individual(s)
    ind$set_state("height", g$height)
    ind$set_state("log_area_leaf_departure", g$phi)
    ind$set_state("log_area_sapwood_departure", g$psi)

    ## Capacity depends on the state, so read it off a seeded twin.
    ref <- TF24_Individual(s)
    ref$set_state("height", g$height)
    ref$set_state("log_area_leaf_departure", g$phi)
    ref$set_state("log_area_sapwood_departure", g$psi)
    ref$set_initial_states(env)
    capacity <- ref$state("storage") / s$pars$a_st3

    ind$set_state("storage", g$storage_frac * capacity)
    ind$compute_rates(env)

    got <- list(dphi = ind$rate("log_area_leaf_departure"),
                dpsi = ind$rate("log_area_sapwood_departure"))
    want <- tf24_departure_rates_r(
      s, g$height, g$phi, g$psi, g$storage_frac * capacity,
      P = ind$aux("net_mass_production_dt"),
      A = ind$aux("competition_effect"))

    label <- sprintf("h=%g phi=%g psi=%g r_frac=%g",
                     g$height, g$phi, g$psi, g$storage_frac)
    expect_equal(got$dphi, want$dphi, tolerance = 1e-10, info = label)
    expect_equal(got$dpsi, want$dpsi, tolerance = 1e-10, info = label)
  }
})
