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
    g1_TF24 = 7.5,
    jmax_25 = 157.44,
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

  # Top-level strategy fields: the pars sub-object plus infrastructure.
  expected_top <- c("pars", "control", "collect_all_auxiliary",
                    "birth_rate_x", "birth_rate_y", "is_variable_birth_rate")

  s <- TF24_Strategy()
  expect_inherits(s, "TF24_Strategy")

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

test_that("TF24 collect_all_auxiliary option", {

  s <- TF24_Strategy()
  p <- TF24_Individual(s)
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

  s <- TF24_Strategy(collect_all_auxiliary=TRUE)
  expect_true(s$collect_all_auxiliary)
  p <- TF24_Individual(s)
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
  s <- TF24_Strategy()
  p <- TF24_Individual(s)

  expect_identical(p$strategy, s)

  ## Set the height to something (here 10)
  h0 <- 10
  p$set_state("height", h0)


  expect_identical(p$state("height"), h0)

  ## Check: Is this redundant now
  ## We now use 
  vars <- p$internals
  expect_identical(p$state("height"), vars$states[which(p$ode_names == "height")])
})



test_that("Critical Names", {
  s <- TF24_Strategy()
  my_names <- TF24_Individual(s)$ode_names
  expect_identical(my_names[1:3], c("height", "mortality", "fecundity"))
})
test_that("TF24_Strategy hyper-parameterisation", {
  s <- TF24_Strategy()

  # lma
  lma <- c(0.1,1)
  ret <- TF24_hyperpar(trait_matrix(lma, "lma"), s)

  expect_true(all(c("lma", "k_l", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "lma"], lma)
  expect_equal(ret[, "k_l"], c(1.46678,0.028600), tolerance=1e-5)
  expect_equal(ret[, "r_l"], c(505.331, 220.633), tolerance=1e-5)

  ## This happens on Linux (and therefore on travis) due to numerical
  ## differences in the integration.
  if ("a_p1" %in% colnames(ret)) {
    a_p1 <- ret[, "a_p1"]
    expect_equal(length(unique(a_p1)), 1L)
    expect_equal(a_p1[[1]], s$pars$a_p1, tolerance=1e-7)
  }

  # wood density
  rho <- c(200,300)
  tf24_hyperpar_rho <- make_TF24_hyperpar(B_hks2 = 1)
  ret <- tf24_hyperpar_rho(trait_matrix(rho, "rho"), s)
  expect_true(all(c("rho", "g1_TF24", "r_s", "r_b") %in% colnames(ret)))
  expect_equal(ret[, "rho"], rho)
  expect_equal(ret[, "g1_TF24"], 7.5 * (rho / 608)^(-1), tolerance = 1e-8)
  expect_equal(ret[, "r_s"], c(20.06000,13.37333), tolerance=1e-5)
  expect_equal(ret[, "r_b"], 2*ret[, "r_s"])

  ## This happens on Linux (and therefore on travis) due to numerical
  ## differences in the integration.
  if ("a_p1" %in% colnames(ret)) {
    a_p1 <- ret[, "a_p1"]
    expect_equal(length(unique(a_p1)), 1L)
    expect_equal(a_p1[[1]], s$pars$a_p1, tolerance=1e-7)
  }

  # vcmax
  vcmax_25 <- c(0, 50,100)
  ret <- TF24_hyperpar(trait_matrix(vcmax_25, "vcmax_25"), s)
  expect_true(all(c("nmass_l", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "vcmax_25"], vcmax_25)
  expect_equal(ret[, "r_l"], c(271.2375, 311.6662, 352.0949), tolerance=1e-5)
  expect_equal(ret[, "nmass_l"], c(0.01234018, 0.01335090, 0.01436162), tolerance=1e-5)
  
  # seed mass
  omega <- 3.8e-5*c(1,2,3)
  ret <- TF24_hyperpar(trait_matrix(omega, "omega"), s)
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
  ret <- TF24_hyperpar(trait_matrix(numeric(0), "lma"), s)
  expect_equal(ret, trait_matrix(numeric(0), "lma"))
})

test_that("TF24_hyperpar no longer produces a_p1/a_p2 and k_I does not affect output", {
  m <- trait_matrix(c(0.1, 0.2), "lma")
  s <- TF24_Strategy()
  ret <- TF24_hyperpar(m, s)

  ## Legacy assimilation parameters a_p1/a_p2 are no longer derived by the
  ## hyperpar function (they remain strategy-level constants, not traits).
  expect_false("a_p1" %in% colnames(ret))
  expect_false("a_p2" %in% colnames(ret))

  ## k_I is no longer used inside the hyperpar function, so varying it in the
  ## strategy does not change the hyperpar output.
  s2 <- TF24_Strategy()
  s2$pars$k_I <- 0.8
  ret2 <- TF24_hyperpar(m, s2)
  expect_equal(ret, ret2)
})

test_that("narea calculation", {
  x <- c(1.38, 3.07, 2.94)
  p0 <- TF24_Parameters()
  m <- trait_matrix(x, "hmat")
  expect_silent(sl <- plant:::generate_strategy(p0, m, hyperpar = TF24_hyperpar, birth_rate = 1.0))

  cmp <- lapply(x, function(xi) generate_strategy(p0, trait_matrix(xi, "hmat"), hyperpar = TF24_hyperpar, birth_rate = 1.0)[[1]])
  expect_equal(sl, cmp)
})

# integration test - runs a full patch meta-population
# the offspring arrival produced integrates all demographic behaviours

test_that("offspring arrival", {

  # This drives a full patch through every demographic process and pins the
  # resulting offspring production so the numbers cannot drift silently.
  #
  # We deliberately use a short patch (max_patch_lifetime = 5) and a low
  # height-at-maturity (hmat = 5) rather than the model defaults. Lowering hmat
  # lets plants mature and reproduce *within* the short patch, so offspring
  # production sits at an O(100) magnitude. That matters for the regression:
  # at the default hmat over so short a patch, reproduction underflows towards
  # zero, and an `expect_equal(tolerance = ...)` against a near-zero target
  # degenerates into a vacuous absolute comparison that any small number passes.
  # The short, low-hmat configuration also cuts TF24's integration time by
  # roughly 7x (TF24 is costly per step -- each rate evaluation runs the
  # root-collar-psi leaf optimisation) while keeping the test a genuine
  # end-to-end check.
  #
  # Tolerance is 1e-3 (relative), not machine-epsilon: offspring_production is an
  # integrated SCM output, so per-step FMA/rounding differences accumulate and it
  # is only reproducible to ~1e-5 across architectures (cf. the arm64 FMA fix in
  # #524). 1e-3 absorbs that while still catching any real regression, which
  # moves these values by whole units or orders of magnitude.
  p0 <- scm_base_parameters("TF24")
  env <- Environment("TF24")
  ctrl <- Control()
  p0$max_patch_lifetime <- 5

  # one species
  p1 <- add_strategies(p0, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                       hyperpar = TF24_hyperpar, birth_rate = list(20))

  out <- run_scm(p1, env, ctrl)
  expect_equal(out$offspring_production, 514.8, tolerance = 1e-3)

  # two species: the second strategy has a moderately higher lma (0.10 vs
  # 0.0825), so it grows more slowly and is partly shaded, but still matures and
  # reproduces within the patch. Both species therefore carry non-zero offspring
  # production -- a genuine two-species coexistence check rather than one
  # strategy collapsing to ~0 (which would make the second value a vacuous
  # near-zero comparison).
  p2 <- add_strategies(p0, trait_matrix(c(0.0825, 0.10, 5, 5), c("lma", "hmat")),
                       hyperpar = TF24_hyperpar, birth_rate = list(20, 20))

  out <- run_scm(p2, env, ctrl)
  expect_equal(out$offspring_production, c(364.02926, 78.28593), tolerance = 1e-3)
})

# Water mass-balance: transpiration integrated up the stem side of every
# individual must match the water depleted from the soil on the root side, under
# a time-varying rainfall driver. Reduced to 5 soil depths (from 15) purely for
# speed (~5x faster, no material change to the closure). The check is
# deliberately one-sided (1 - ratio < tol, i.e. ratio > 1 - tol): over so short
# a transient patch the cumulative-flux closure does not settle to a tight
# two-sided tolerance. See #533 for tightening it.

test_that("E conservation", {

max_patch_lifetime <-2
p0 <- scm_base_parameters("TF24", "TF24_Env")
p0$max_patch_lifetime <- max_patch_lifetime
traits <- trait_matrix(c(0.07), c("lma"))
p1 <- add_strategies(p0, traits)

env <- Environment("TF24")
env$set_soil_number_of_depths(5)
env$set_soil_water_state(rep(c(0.2), times = 5))
x = seq(0,max_patch_lifetime,length.out = 100)
y = 0.25*sin(2*pi*x) + 1
env$extrinsic_drivers_set_variable("rainfall", x=x, y=y)
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

expect_true(1 - (stem_side/root_side$root_side[-1])[length(stem_side)] < 5e-2)
})

test_that("SCM cohort-density blow-up fails gracefully (#550)", {
  # Extreme seasonal drought drives the SCM size-density (characteristic)
  # equations to run away: a cohort density overflows to +Inf, or the
  # density-weighted resource uptake drives a soil-water state non-finite.
  # Either way run_scm must abort with an actionable message naming the SCM
  # size-density equations / environment state, not an opaque downstream error
  # ("Detected non-finite contribution", "non-finite psi_soil"). This blows up
  # early in the run, so it is cheap despite being a full TF24 run.
  mpl <- 30
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- mpl
  p1 <- add_strategies(p0, trait_matrix(0.07, "lma"))

  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.2, 5))
  x <- seq(0, mpl, length.out = mpl * 6)
  y <- 0.4 * sin(2 * pi * x) + 0.5   # rainfall sweeps [0.1, 0.9]
  env$extrinsic_drivers_set_variable("rainfall", x = x, y = y)

  expect_error(run_scm(p1, env),
               "SCM size-density|Non-finite environment state")
})

