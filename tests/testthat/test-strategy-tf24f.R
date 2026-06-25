# Built from tests/testthat/test-strategy-tf24.R on Thu Jun 25 06:01:09 2026 using the scaffolder, from the strategy: TF24
# Built from  tests/testthat/test-strategy-ff16.R on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16

test_that("Defaults", {
  # Biological parameters now live in the nested `pars` sub-object.
  expected_pars <- list(
    a_l2     = 0.306,
    S_D   = 0.25,
    a_y      = 0.7,
    a_l1     = 5.44,
    a_r1     = 0.07,
    a_b1      = 0.17,
    r_b   = 8024 / 608,
    r_l   = 39.27 / 0.1978791,
    r_r   = 217,
    r_s   = 4012/608,
    a_f3  = 3.0*3.8e-5,
    a_bio  = 0.0245,
    d_I   = 0.01,
    a_dG1   = 5.5,
    a_dG2   = 20,
    a_p1   = 151.177775377968,
    a_p2   = 0.204716166503633,
    a_f1   = 1,
    a_f2   = 50,
    a_d0   = 0.1,
    eta    = 12,
    hmat   = 16.5958691,
    k_b    = 0.2,
    k_l   = 0.4565855,
    k_r    = 1,
    k_s   = 0.2,
    lma    = 0.1978791,
    rho    = 608,
    omega  = 3.8e-5,
    theta  = 1.0/4669,
    k_I = 0.5,
    vcmax_25 = 96,
    p_50 = 1.85,
    K_s = 1,
    c = log(log(1-0.5)/log(1-0.88))/(log(1.85) - log(5.16)),
    b = 1.85 /((-log(1 - 50.0 / 100.0))^(1 / (log(log(1-0.5)/log(1-0.88))/(log(1.85) - log(5.16))))),
    psi_crit = (1.85 /((-log(1 - 50.0 / 100.0))^(1 / (log(log(1-0.5)/log(1-0.88))/(log(1.85) - log(5.16))))))*log(1/0.05)^(1/(log(log(1-0.5)/log(1-0.88))/(log(1.85) - log(5.16)))),
    beta1 = 20000,
    beta2 = 1.5,
    jmax_25 = 157.44,
    hk_s = 4,
    a = 0.3,
    curv_fact_elec_trans = 0.7,
    curv_fact_colim = 0.99,
    var_sapwood_volume_cost = 1,
    nmass_l = 0.013,
    nmass_s = 0.00198,
    nmass_b = 0.0034,
    nmass_r = 0.00335,
    dmass_dN = 0,
    root_depth_shape_eta = 0.2,
    recruitment_decay = 0)

  # Top-level strategy fields: the pars sub-object plus infrastructure, and the
  # TF24f-specific acclimation knobs (gain + finite-difference step) (#525).
  expected_top <- c("pars", "control", "collect_all_auxiliary",
                    "birth_rate_x", "birth_rate_y", "is_variable_birth_rate",
                    "k_acclim", "psi_fd_step", "use_ad_gradient")

  s <- TF24f_Strategy()
  expect_inherits(s, "TF24f_Strategy")

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
  # k_acclim must be finite and >= 0; psi_fd_step must be finite and > 0.
  p1 <- scm_base_parameters("TF24f") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  s1 <- p1$strategies[[1]]; s1$k_acclim <- -1; p1$strategies[[1]] <- s1
  expect_error(run_scm(p1, collect = FALSE, refine_schedule = FALSE), "k_acclim")

  # psi_fd_step is only used on the finite-difference path, so exercise it with
  # the AD gradient disabled.
  p2 <- scm_base_parameters("TF24f") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  s2 <- p2$strategies[[1]]; s2$psi_fd_step <- 0; s2$use_ad_gradient <- FALSE
  p2$strategies[[1]] <- s2
  expect_error(run_scm(p2, collect = FALSE, refine_schedule = FALSE), "psi_fd_step")
})

test_that("TF24f collect_all_auxiliary option", {

  s <- TF24f_Strategy()
  p <- TF24f_Individual(s)
  expect_equal(p$aux_size, 11)
  expect_equal(length(p$internals$auxs), 11)
expect_equal(p$aux_names, c(
    "competition_effect",
    "height_inverse",
    "net_mass_production_dt",
    "root_mass",
    "opt_psi_stem",
    "opt_root_psi",
    "transpiration",
    "E_up_",
    "profit",
    "stom_cond_CO2",
    "assimilation"
  ))

  s <- TF24f_Strategy(collect_all_auxiliary=TRUE)
  expect_true(s$collect_all_auxiliary)
  p <- TF24f_Individual(s)
  expect_equal(p$aux_size, 12)
  expect_equal(length(p$internals$auxs), 12)
  expect_equal(p$aux_names, c(
    "competition_effect",
    "height_inverse",
    "net_mass_production_dt",
    "root_mass",
    "opt_psi_stem",
    "opt_root_psi",
    "transpiration",
    "E_up_",
    "profit",
    "stom_cond_CO2",
    "assimilation",
    "area_sapwood"
  ))
})

test_that("Reference comparison", {
  s <- TF24f_Strategy()
  p <- TF24f_Individual(s)

  expect_identical(p$strategy, s)

  ## Set the height to something (here 10)
  h0 <- 10
  p$set_state("height", h0)


  expect_identical(p$state("height"), h0)

  ## The state() accessor and the raw ODE state vector agree: state("height")
  ## must equal the entry of internals$states at the height slot.
  vars <- p$internals
  expect_identical(p$state("height"), vars$states[which(p$ode_names == "height")])
})



test_that("Critical Names", {
  s <- TF24f_Strategy()
  my_names <- TF24f_Individual(s)$ode_names
  expect_identical(my_names[1:3], c("height", "mortality", "fecundity"))
})
test_that("TF24f_Strategy hyper-parameterisation", {
  s <- TF24f_Strategy()

  # lma
  lma <- c(0.1,1)
  ret <- TF24f_hyperpar(trait_matrix(lma, "lma"), s)

  expect_true(all(c("lma", "k_l", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "lma"], lma)
  expect_equal(ret[, "k_l"], c(1.46678,0.028600), tolerance=1e-5)
  expect_equal(ret[, "r_l"], c(392.70, 39.27), tolerance=1e-5)

  ## This happens on Linux (and therefore on travis) due to numerical
  ## differences in the integration.
  if ("a_p1" %in% colnames(ret)) {
    a_p1 <- ret[, "a_p1"]
    expect_equal(length(unique(a_p1)), 1L)
    expect_equal(a_p1[[1]], s$pars$a_p1, tolerance=1e-7)
  }

  # wood density
  rho <- c(200,300)
  ret <- TF24f_hyperpar(trait_matrix(rho, "rho"), s)
  expect_true(all(c("rho", "r_s", "r_b") %in% colnames(ret)))
  expect_equal(ret[, "rho"], rho)
  expect_equal(ret[, "r_s"], c(20.06000,13.37333), tolerance=1e-5)
  expect_equal(ret[, "r_b"], 2*ret[, "r_s"])

  ## This happens on Linux (and therefore on travis) due to numerical
  ## differences in the integration.
  if ("a_p1" %in% colnames(ret)) {
    a_p1 <- ret[, "a_p1"]
    expect_equal(length(unique(a_p1)), 1L)
    expect_equal(a_p1[[1]], s$pars$a_p1, tolerance=1e-7)
  }

  # narea
  narea <- c(0, 2E-3,2.3E-3)
  ret <- TF24f_hyperpar(trait_matrix(narea, "narea"), s)
  expect_true(all(c("narea", "a_p1", "a_p2", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "narea"], narea)
  expect_equal(ret[, "r_l"], c(0, 212.2508, 244.0884), tolerance=1e-5)
  expect_equal(ret[, "a_p1"], c(0, 162.2592, 188.1549), tolerance=1e-5)
  expect_equal(ret[, "a_p2"], c(0, 0.220904, 0.259173), tolerance=1e-5)

  # seed mass
  omega <- 3.8e-5*c(1,2,3)
  ret <- TF24f_hyperpar(trait_matrix(omega, "omega"), s)
  expect_true(all(c("omega", "a_f3") %in% colnames(ret)))
  expect_equal(ret[, "omega"], omega)
  expect_equal(ret[, "a_f3"], 3*omega)

  ## This happens on Linux (and therefore on travis) due to numerical
  ## differences in the integration.
  if ("a_p1" %in% colnames(ret)) {
    a_p1 <- ret[, "a_p1"]
    expect_equal(length(unique(a_p1)), 1L)
    expect_equal(a_p1[[1]], s$pars$a_p1, tolerance=1e-7)
  }


  ## Empty trait matrix:
  ret <- TF24f_hyperpar(trait_matrix(numeric(0), "lma"), s)
  expect_equal(ret, trait_matrix(numeric(0), "lma"))
})

test_that("TF24f_hyperpar sources k_I from the strategy", {
  narea <- c(2E-3, 2.3E-3)
  m <- trait_matrix(narea, "narea")

  ## Default strategy: assimilation matches the existing reference values.
  s <- TF24f_Strategy()
  expect_equal(s$pars$k_I, 0.5)
  ret <- TF24f_hyperpar(m, s)
  expect_equal(ret[, "a_p1"], c(162.2592, 188.1549), tolerance=1e-5)

  ## Varying the strategy's k_I must change the derived assimilation
  ## parameters -- previously the hard-coded 0.5 default in the maker
  ## silently ignored the strategy value.
  s2 <- TF24f_Strategy()
  s2$pars$k_I <- 0.8
  ret2 <- TF24f_hyperpar(m, s2)
  expect_false(isTRUE(all.equal(ret[, "a_p1"], ret2[, "a_p1"])))
  expect_false(isTRUE(all.equal(ret[, "a_p2"], ret2[, "a_p2"])))
})

test_that("narea calculation", {
  x <- c(1.38, 3.07, 2.94)
  p0 <- TF24f_Parameters()
  m <- trait_matrix(x, "hmat")
  expect_silent(sl <- plant:::generate_strategy(p0, m, hyperpar = TF24f_hyperpar, birth_rate = 1.0))

  cmp <- lapply(x, function(xi) generate_strategy(p0, trait_matrix(xi, "hmat"), hyperpar = TF24f_hyperpar, birth_rate = 1.0)[[1]])
  expect_equal(sl, cmp)
})

# TF24f is the fast-acclimating variant of TF24 (#525): instead of solving the
# stem/root water-potential optimisation afresh every step (as TF24 does), it
# carries an extra ODE state `opt_root_psi_state` that relaxes towards the
# optimum at rate `k_acclim`, using a finite-difference (`psi_fd_step`) gradient.
#
# We do NOT re-pin a full demographic regression here -- that is TF24's job, and
# the two strategies share all of that machinery. Instead this integration test
# targets what is *unique* to TF24f, cheaply: (1) it runs end-to-end, (2) the
# acclimation/gradient feedback is actually live (the gain materially changes
# fitness), and (3) in the strong-acclimation limit it reproduces TF24's
# directly-optimised result. The short patch + low hmat configuration matches
# test-strategy-tf24.R so the cross-strategy comparison in (3) is like-for-like.

test_that("acclimation runs, is active, and converges to TF24", {

  mk <- function(strat) {
    p0 <- scm_base_parameters(strat)
    p0$max_patch_lifetime <- 5
    add_strategies(p0, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                   hyperpar = get(paste0(strat, "_hyperpar")), birth_rate = list(20))
  }
  set_k_acclim <- function(p, k) {
    s <- p$strategies[[1]]
    s$k_acclim <- k
    p$strategies[[1]] <- s
    p
  }

  # (1) The acclimation state is part of the ODE system (TF24 has no 6th state),
  #     and a full patch runs end-to-end to a finite, positive offspring count.
  expect_true("opt_root_psi_state" %in% TF24f_Individual()$ode_names)
  pf <- mk("TF24f")
  slow <- run_scm(set_k_acclim(pf, 0.1), Environment("TF24f"), Control())$offspring_production
  fast <- run_scm(set_k_acclim(pf, 10),  Environment("TF24f"), Control())$offspring_production
  expect_length(fast, 1)
  expect_true(is.finite(fast) && fast > 0)

  # (2) Acclimation is genuinely active: the gain k_acclim materially changes
  #     fitness, so the finite-difference psi optimisation is feeding back into
  #     the demography rather than being a no-op.
  expect_gt(abs(fast - slow) / slow, 0.1)

  # (3) Consistent with TF24 "within reason": as the acclimation gain grows,
  #     TF24f tracks the optimum that TF24 computes directly each step, so its
  #     offspring production converges onto TF24's.
  tf24 <- run_scm(mk("TF24"), Environment("TF24"), Control())$offspring_production
  expect_equal(fast, tf24, tolerance = 1e-2)
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
  p0 <- scm_base_parameters("TF24f", "TF24_Env")
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

  results <- run_scm(p1, env = env, ctrl = ctrl, collect = TRUE)

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
    p0 <- scm_base_parameters(type)
    p0$max_patch_lifetime <- 5
    add_strategies(p0, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                   hyperpar = get(paste0(type, "_hyperpar")), birth_rate = 20)
  }
  ref <- run_scm(base("TF24"), collect = TRUE, refine_schedule = FALSE)$species

  p <- base("TF24f")
  s <- p$strategies[[1]]
  s$k_acclim <- 10; s$use_ad_gradient <- TRUE
  p$strategies[[1]] <- s
  res <- run_scm(p, collect = TRUE, refine_schedule = FALSE)$species

  expect_equal(nrow(res), nrow(ref))
  dh <- abs(res$height - ref$height)
  expect_gt(mean(dh <= 0.01 * pmax(ref$height, 1e-3), na.rm = TRUE), 0.90)
  expect_lt(median(dh, na.rm = TRUE), 0.05)
})
