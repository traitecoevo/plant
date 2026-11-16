test_that("Defaults", {
  # TF24f carries TF24's parameter struct, so its defaults are TF24's, which
  # test-model-tf24.R pins value by value.
  expected_pars <- unclass(TF24_Model()$pars)

  # Top-level strategy fields: the pars sub-object plus infrastructure, and
  # TF24f's acclimation gain.
  expected_top <- c("pars", "control", "collect_all_auxiliary",
                    "birth_rate_x", "birth_rate_y", "is_variable_birth_rate",
                    "k_acclim")

  s <- TF24f_Model()
  expect_inherits(s, "TF24f_Model")

  expect_identical(sort(names(s)), sort(expected_top))
  expect_identical(s$control, Control())
  expect_identical(s$collect_all_auxiliary, FALSE)
  expect_identical(s$birth_rate_x, numeric(0))
  expect_identical(s$birth_rate_y, c(1.0))
  expect_identical(s$is_variable_birth_rate, FALSE)

  pars_keys <- sort(names(expected_pars))
  expect_identical(sort(names(s$pars)), pars_keys)
  expect_identical(unclass(s$pars)[pars_keys], expected_pars[pars_keys])
})

test_that("TF24f validates acclimation knobs", {
  # k_acclim must be finite and >= 0.
  p1 <- testing_parameters("TF24f") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  s1 <- p1$strategies[[1]]; s1$k_acclim <- -1; p1$strategies[[1]] <- s1
  expect_error(run(p1), "k_acclim")
})

test_that("TF24f collect_all_auxiliary option", {

  s <- TF24f_Model()
  p <- TF24f_Individual(s)
  expect_equal(p$aux_size, 15)
  expect_equal(length(p$internals$auxs), 15)
expect_equal(p$aux_names, c(
    "competition_effect",
    "height_inverse",
    "net_mass_production_dt",
    "root_mass",
    "opt_psi_leaf",
    "opt_root_psi",
    "transpiration",
    "E_up_",
    "profit",
    "shadow_cost",
    "stom_cond_CO2",
    "assimilation",
    "Tleaf", "leaf_marginal_return", "sapwood_marginal_return"
  ))

  s <- TF24f_Model(collect_all_auxiliary=TRUE)
  expect_true(s$collect_all_auxiliary)
  p <- TF24f_Individual(s)
  expect_equal(p$aux_size, 16)
  expect_equal(length(p$internals$auxs), 16)
  expect_equal(p$aux_names, c(
    "competition_effect",
    "height_inverse",
    "net_mass_production_dt",
    "root_mass",
    "opt_psi_leaf",
    "opt_root_psi",
    "transpiration",
    "E_up_",
    "profit",
    "shadow_cost",
    "stom_cond_CO2",
    "assimilation",
    "Tleaf", "leaf_marginal_return", "sapwood_marginal_return",
    "area_sapwood"
  ))
})

test_that("TF24f_Model hyper-parameterisation", {
  s <- TF24f_Model()

  # lma
  lma <- c(0.1,1)
  ret <- default_hyperpar("TF24")(trait_matrix(lma, "lma"), s)

  expect_true(all(c("lma", "k_l", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "lma"], lma)
  expect_equal(ret[, "k_l"], c(1.46678,0.028600), tolerance=1e-5)
  expect_equal(ret[, "r_l"], c(505.33080, 220.63308), tolerance=1e-5)

  # wood density
  rho <- c(200,300)
  tf24f_hyperpar_rho <- default_hyperpar("TF24", B_hks2 = 1)
  ret <- tf24f_hyperpar_rho(trait_matrix(rho, "rho"), s)
  expect_true(all(c("rho", "TF24_cost_scale", "r_s", "r_b") %in% colnames(ret)))
  expect_equal(ret[, "rho"], rho)
  expect_equal(ret[, "TF24_cost_scale"], 7.5 * (rho / 608)^(-1), tolerance = 1e-8)
  expect_equal(ret[, "r_s"], c(20.06000,13.37333), tolerance=1e-5)
  expect_equal(ret[, "r_b"], 2*ret[, "r_s"])

  # sapwood conductivity and hydraulic traits
  K_s <- c(1, 2, 3)
  ret <- default_hyperpar("TF24")(trait_matrix(K_s, "K_s"), s)
  expect_true(all(c("K_s", "stem_P50", "stem_b", "psi_crit") %in% colnames(ret)))
  # B_c2 = 0 makes stem_c the default at every K_s, so the filter drops it.
  expect_false("stem_c" %in% colnames(ret))
  expect_equal(ret[, "K_s"], K_s)
  expect_true(all(diff(ret[, "stem_P50"]) < 0))

  # vcmax_25 contributes to leaf N and therefore respiration
  vcmax_25 <- c(0, 96, 160)
  ret <- default_hyperpar("TF24")(trait_matrix(vcmax_25, "vcmax_25"), s)
  expect_true(all(c("vcmax_25", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "vcmax_25"], vcmax_25)
  expect_true(all(diff(ret[, "r_l"]) > 0))

  # seed mass
  omega <- 3.8e-5*c(1,2,3)
  ret <- default_hyperpar("TF24")(trait_matrix(omega, "omega"), s)
  expect_true(all(c("omega", "a_f3") %in% colnames(ret)))
  expect_equal(ret[, "omega"], omega)
  expect_equal(ret[, "a_f3"], 3*omega)

  expect_false(any(c("a_p1", "a_p2") %in% colnames(ret)))


  ## Empty trait matrix:
  ret <- default_hyperpar("TF24")(trait_matrix(numeric(0), "lma"), s)
  expect_equal(ret, trait_matrix(numeric(0), "lma"))
})

test_that("the TF24 trade-off function sources k_I from the strategy", {
  m <- trait_matrix(c(96, 120), "vcmax_25")

  s <- TF24f_Model()
  expect_equal(s$pars$k_I, 0.5)
  ret <- default_hyperpar("TF24")(m, s)

  s2 <- TF24f_Model()
  s2$pars$k_I <- 0.8
  ret2 <- default_hyperpar("TF24")(m, s2)

  # TF24f now follows TF24's updated hyperparameterisation, which no longer
  # computes the legacy assimilation-fit columns a_p1/a_p2 from k_I.
  expect_equal(ret, ret2, tolerance = 1e-12)
})

test_that("narea calculation", {
  x <- c(1.38, 3.07, 2.94)
  p0 <- TF24f_Parameters()
  m <- trait_matrix(x, "hmat")
  expect_silent(sl <- plant:::generate_strategy(p0, m, hyperpar = default_hyperpar("TF24"), birth_rate = 1.0))

  cmp <- lapply(x, function(xi) generate_strategy(p0, trait_matrix(xi, "hmat"), hyperpar = default_hyperpar("TF24"), birth_rate = 1.0)[[1]])
  expect_equal(sl, cmp)
})

# TF24f is TF24 with a plastic hydraulic operating point: instead of solving the
# leaf's optimal root-collar potential afresh at every step (as TF24 does), it
# carries an extra ODE state `collar_position`, the logit of the collar's
# position between the wet and dry ends of its feasible interval, which climbs
# the exact gradient of profit at the gain `k_acclim`.
#
# We do NOT re-pin a full demographic regression here -- that is TF24's job, and
# the two strategies share all of that machinery. Instead this integration test
# targets what is *unique* to TF24f, cheaply: (1) it runs end-to-end, (2) the
# acclimation/gradient feedback is actually live (the gain materially changes
# fitness), and (3) in the strong-acclimation limit it reproduces TF24's
# directly-optimised result. The short patch + low hmat configuration matches
# test-model-tf24.R so the cross-strategy comparison in (3) is like-for-like.

tf24f_stand <- function(strat, lifetime = 5, lma = 0.0825, k_acclim = NULL) {
  p <- testing_parameters(strat)
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(c(lma, 5), c("lma", "hmat")),
                      hyperpar = default_hyperpar(strat), birth_rate = list(20))
  if (!is.null(k_acclim)) {
    s <- p$strategies[[1]]
    s$k_acclim <- k_acclim
    p$strategies[[1]] <- s
  }
  p
}

# A five-layer soil under a rain that moves the feasible interval through the
# year: dry pulses, where profit still rises at the interval's dry end for half
# of each year, or a seasonal rain, where shaded plants sit at its wet end.
tf24f_soil <- function(strat, rain, lifetime = 10) {
  env <- Environment(strat)
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.2, 5))
  x <- seq(0, lifetime, length.out = if (rain == "dry pulses") 400 else 100)
  y <- if (rain == "dry pulses") 0.3 * pmax(sin(2 * pi * x), 0)
       else 0.25 * sin(2 * pi * x) + 1
  env$extrinsic_drivers_set_variable("rainfall", x = x, y = y)
  env
}

test_that("acclimation runs, is active, and converges to TF24", {

  # (1) The acclimation state is part of the ODE system (TF24 has no such
  #     state), and a full patch runs end-to-end to a finite, positive
  #     offspring count.
  expect_true("collar_position" %in% TF24f_Individual()$ode_names)
  at_gain <- function(k) {
    run(tf24f_stand("TF24f", k_acclim = k), Environment("TF24f"),
        Control())$offspring_production
  }
  slow <- at_gain(0.01)
  fast <- at_gain(10)
  expect_length(fast, 1)
  expect_true(is.finite(fast) && fast > 0)

  # (2) Acclimation is genuinely active: the gain k_acclim materially changes
  #     fitness, so the tracked position is feeding back into the demography
  #     rather than being a no-op.
  #
  #     ⚠️ Offspring production is NOT monotone in k_acclim, so a pair of gains
  #     can sit at the same height on either side of a turning point and differ
  #     by little with the feedback fully live. Measured here at hmat = 5,
  #     max_patch_lifetime = 5, on the birth-date coordinate TF24f defaults to:
  #
  #       k_acclim    0.01     0.1      1       10      100
  #       offspring   367.4   314.4   320.4   330.0   330.7
  #
  #     A plant that barely adjusts keeps the position it was born with, and
  #     differs from one that tracks by 11%.
  expect_gt(abs(fast - slow) / slow, 0.05)

  # (3) As the acclimation gain grows, TF24f tracks the optimum that TF24
  #     computes directly each step, so its offspring production converges
  #     onto TF24's. Measured gap against TF24 as the gain rises:
  #
  #       k_acclim     1       3       10      30      100
  #       rel. gap   3.1%    0.98%   0.26%   0.083%  0.024%
  #
  #     TF24 is integrated to 1e-6 here. Monotone to zero from k = 1, which is
  #     what makes this a tracking lag and not a structural divergence between
  #     the two strategies.
  tight <- Control()
  tight$ode_tol_rel <- 1e-6
  tight$ode_tol_abs <- 1e-6
  tf24 <- run(tf24f_stand("TF24"), Environment("TF24"), tight)$offspring_production
  converged <- at_gain(100)
  expect_equal(converged, tf24, tolerance = 1e-3)
  # ...and the approach is monotone: a higher gain is never further away.
  expect_lt(abs(converged - tf24), abs(fast - tf24))
})

test_that("under drought the collar stays in its interval and follows TF24", {
  # Under dry pulses profit still rises at the dry end of the feasible interval
  # for about half of each year, and the interval itself moves by megapascals
  # as the soil dries and rewets. The position is the collar's place between
  # the ends, so it moves with them, and the barrier holds it a short step
  # inside the dry end: a finite position throughout, near logit(1 - 1e-4),
  # where no state held in MPa stayed inside the interval. Measured gap in
  # offspring production against TF24 (27.38):
  #
  #   k_acclim     1      10      30     100
  #   rel. gap   4.0%    4.9%    1.6%   0.45%
  lifetime <- 10
  at_gain <- function(k) {
    run(tf24f_stand("TF24f", lifetime, lma = 0.07, k_acclim = k),
        tf24f_soil("TF24f", "dry pulses"), Control(), record = "steps")
  }
  tf24 <- run(tf24f_stand("TF24", lifetime, lma = 0.07),
              tf24f_soil("TF24", "dry pulses"), Control())$offspring_production
  slow <- at_gain(1)
  position <- tidy(slow)$species$collar_position
  position <- position[!is.na(position)]
  expect_true(all(is.finite(position)))
  expect_gt(max(position), 8)
  expect_lt(max(position), 12)
  expect_equal(slow$offspring_production, tf24, tolerance = 0.08)
  expect_equal(at_gain(30)$offspring_production, tf24, tolerance = 0.03)
})

test_that("offspring production is smooth in a trait at the interval's ends", {
  # Under a seasonal rain shaded plants sit at the wet end of the interval,
  # where no water flows. A state clamped to that end read the no-flow marginal,
  # an exact zero beside a finite value, so offspring production jumped as a
  # trait moved: its response to lma, differenced over runs pinned to the
  # recorded steps, went from -1566 to -301,097 between relative steps of 1e-7
  # and 1e-6. A smooth model gives the same response at every step.
  lifetime <- 10
  ctrl <- Control()
  env <- function() tf24f_soil("TF24f", "seasonal")
  base <- run(tf24f_stand("TF24f", lifetime, lma = 0.07), env(), ctrl,
              record = "states")$solver
  pinned <- base$parameters
  pinned$ode_times <- base$ode_times
  pinned$ode_step_sizes <- base$ode_step_sizes
  pinned$ode_soil_input_slopes <- base$ode_soil_input_slopes
  pinned$ode_soil_substep_ends <- base$ode_soil_substep_ends
  response <- function(rel) {
    q <- pinned
    q$strategies[[1]]$pars$lma <- q$strategies[[1]]$pars$lma * (1 + rel)
    (run(q, env(), ctrl)$offspring_production -
       base$patch$offspring_production) / rel
  }
  expect_equal(response(1e-6), response(1e-4), tolerance = 1e-3)
})

# Water mass-balance: the transpiration integrated up the stem side of every
# individual must match the water depleted from the soil on the root side. This
# exercises the coupled soil-water / hydraulics machinery end-to-end under a
# time-varying rainfall driver.
#
# Reduced to 5 soil depths (from 15) purely for speed -- ~5x faster with no
# material change to the closure.
#
# Two-sided check (#533): the original assertion was one-sided (1 - ratio < tol),
# which could only catch the stem side falling *below* the root side, not an
# overshoot. The overshoot it would otherwise have masked is a short-patch
# transient, not a conservation gap: the end-of-patch closure converges to 1.0 as
# the patch lengthens (ratio ~1.25 at lifetime 2 -> ~1.01 at lifetime 10 -> ~1.00
# at lifetime 40), because the size distribution and the tracked leaf states are
# still equilibrating early in a short patch. We therefore run a longer patch
# (lifetime 10, ~2.4s) where instantaneous flux balance has settled and assert a
# genuine two-sided closure |1 - ratio| < tol.

test_that("E conservation", {

  max_patch_lifetime <- 10
  p0 <- testing_parameters("TF24f")
  p0$max_patch_lifetime <- max_patch_lifetime
  traits <- trait_matrix(c(0.07), c("lma"))
  p1 <- add_strategies(p0, traits)

  env <- Environment("TF24f")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(c(0.2), times = 5))
  x <- seq(0, max_patch_lifetime, length.out = 100)
  y <- 0.25 * sin(2 * pi * x) + 1
  env$extrinsic_drivers_set_variable("rainfall", x = x, y = y)
  ctrl <- Control()

  results <- tidy(run(p1, env = env, control = ctrl, record = "steps"))

  results %>%
    expand_state() %>%
    purrr::pluck("species") %>%
    dplyr::mutate(E_indiv = E_up_ * area_leaf * 60 * 60 * 12 * 365 / 1000) %>%
    integrate_over_size_distribution() %>%
    dplyr::pull(E_indiv) -> stem_side

  results$env$soil_moist_cumulative_flux %>%
    dplyr::mutate(
      root_side = (sum_resource_depletion - dplyr::lag(sum_resource_depletion)) /
                  (time - dplyr::lag(time))) -> root_side

  closure <- (stem_side / root_side$root_side[-1])[length(stem_side)]
  expect_equal(closure, 1, tolerance = 5e-2)
})

test_that("TF24f AD gradient tracks TF24 (#527)", {
  # End-to-end regression: with the exact AD/IFT gradient, a TF24f patch tracks
  # the TF24 (full-optimisation) patch closely at a moderate acclimation gain.
  #
  # Short patch + low hmat (matching the offspring-arrival / acclimation tests):
  # the tracking property is about per-step optimisation accuracy, not patch
  # length, so a 5-year patch reproduces it (cohorts agree to <1% at every step)
  # while running ~35x faster than the default horizon.
  base <- function(type) {
    p0 <- testing_parameters(type)
    p0$max_patch_lifetime <- 5
    add_strategies(p0, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                   hyperpar = default_hyperpar(type), birth_rate = 20)
  }
  ref <- tidy(run(base("TF24"), record = "steps"))$species

  p <- base("TF24f")
  s <- p$strategies[[1]]
  s$k_acclim <- 10
  p$strategies[[1]] <- s
  res <- tidy(run(p, record = "steps"))$species

  expect_equal(nrow(res), nrow(ref))
  dh <- abs(res$height - ref$height)
  expect_gt(mean(dh <= 0.01 * pmax(ref$height, 1e-3), na.rm = TRUE), 0.90)
  expect_lt(median(dh, na.rm = TRUE), 0.05)
})
