
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
    a_f4   = 0,
    a_f5   = 0,
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
    recruitment_decay = 0)

  # Top-level strategy fields: the pars sub-object plus infrastructure.
  expected_top <- c("pars", "control", "collect_all_auxiliary",
                    "birth_rate_x", "birth_rate_y", "is_variable_birth_rate")

  s <- FF16_Strategy()
  expect_inherits(s, "FF16_Strategy")

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

test_that("FF16 collect_all_auxiliary option", {

  s <- FF16_Strategy()
  p <- FF16_Individual(s)
  expect_equal(p$aux_size, 3)
  expect_equal(length(p$internals$auxs), 3)
  expect_equal(p$aux_names, c(
    "competition_effect",
    "height_inverse",
    "net_mass_production_dt"
  ))

  s <- FF16_Strategy(collect_all_auxiliary=TRUE)
  expect_true(s$collect_all_auxiliary)
  p <- FF16_Individual(s)
  expect_equal(p$aux_size, 4)
  expect_equal(length(p$internals$auxs), 4)
  expect_equal(p$aux_names, c(
    "competition_effect",
    "height_inverse",
    "net_mass_production_dt",
    "area_sapwood"
  ))
})

test_that("Reference comparison", {
  s <- FF16_Strategy()
  p <- FF16_Individual(s)

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
  s <- FF16_Strategy()
  my_names <- FF16_Individual(s)$ode_names
  expect_identical(my_names[1:3], c("height", "mortality", "fecundity"))
})

test_that("FF16_Strategy hyper-parameterisation", {
  s <- FF16_Strategy()

  # lma
  lma <- c(0.1,1)
  ret <- FF16_hyperpar(trait_matrix(lma, "lma"), s)

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
  ret <- FF16_hyperpar(trait_matrix(rho, "rho"), s)
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
  ret <- FF16_hyperpar(trait_matrix(narea, "narea"), s)
  expect_true(all(c("narea", "a_p1", "a_p2", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "narea"], narea)
  expect_equal(ret[, "r_l"], c(0, 212.2508, 244.0884), tolerance=1e-5)
  expect_equal(ret[, "a_p1"], c(0, 162.2592, 188.1549), tolerance=1e-5)
  expect_equal(ret[, "a_p2"], c(0, 0.220904, 0.259173), tolerance=1e-5)

  # seed mass
  omega <- 3.8e-5*c(1,2,3)
  ret <- FF16_hyperpar(trait_matrix(omega, "omega"), s)
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
  ret <- FF16_hyperpar(trait_matrix(numeric(0), "lma"), s)
  expect_equal(ret, trait_matrix(numeric(0), "lma"))
})

test_that("FF16_hyperpar sources k_I from the strategy", {
  narea <- c(2E-3, 2.3E-3)
  m <- trait_matrix(narea, "narea")

  ## Default strategy: assimilation matches the existing reference values.
  s <- FF16_Strategy()
  expect_equal(s$pars$k_I, 0.5)
  ret <- FF16_hyperpar(m, s)
  expect_equal(ret[, "a_p1"], c(162.2592, 188.1549), tolerance=1e-5)

  ## Varying the strategy's k_I must change the derived assimilation
  ## parameters -- previously the hard-coded 0.5 default in the maker
  ## silently ignored the strategy value.
  s2 <- FF16_Strategy()
  s2$pars$k_I <- 0.8
  ret2 <- FF16_hyperpar(m, s2)
  expect_false(isTRUE(all.equal(ret[, "a_p1"], ret2[, "a_p1"])))
  expect_false(isTRUE(all.equal(ret[, "a_p2"], ret2[, "a_p2"])))
})

test_that("narea calculation", {
  x <- c(1.38, 3.07, 2.94)
  p0 <- FF16_Parameters()
  m <- trait_matrix(x, "hmat")
  expect_silent(sl <- plant:::generate_strategy(p0, m, hyperpar = FF16_hyperpar, birth_rate = 1.0))

  cmp <- lapply(x, function(xi) generate_strategy(p0, trait_matrix(xi, "hmat"), hyperpar = FF16_hyperpar, birth_rate = 1.0)[[1]])
  expect_equal(sl, cmp)
})

# integration test - runs a full patch meta-population
# the offspring arrival produced integrates all demographic behaviours
test_that("offspring arrival", {

  p0 <- scm_base_parameters("FF16")
  env <- Environment("FF16")
  ctrl <- Control()
  
  # one species
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))

  out <- run_scm(p1, env, ctrl)
  expect_equal(out$offspring_production, 16.8895016, tolerance=1e-4)
  expect_equal(out$ode_times[c(10, 100)], c(0.000070, 4.215899), tolerance=1e-5)

  # two species
  p2 <- add_strategies(p0, trait_matrix(c(0.0825, 0.2625), "lma"), hyperpar = FF16_hyperpar, birth_rate = list(11.99177, 16.51006))
  
  out <- run_scm(p2, env, ctrl)
  expect_equal(out$offspring_production, c(11.995204, 16.474988), tolerance=1e-5)
  expect_equal(length(out$ode_times), 293)
})

# ---------------------------------------------------------------------------
# Reproductive allocation as a reaction norm on light (a_f4, a_f5)
#
#   RA(H, L) = a_f1 * L^a_f5 / (1 + exp(a_f2 * (1 - H / hmat_eff)))
#   hmat_eff = hmat * (1 + a_f4 * (1 - L))
#
# with L = canopy openness at the top of the plant's own crown.

## Analytic reaction norm, mirroring FF16_Strategy::fraction_allocation_reproduction.
ff16_ra_expected <- function(pars, height, openness) {
  hmat_eff <- pars$hmat * (1 + pars$a_f4 * (1 - openness))
  pars$a_f1 * openness^pars$a_f5 /
    (1 + exp(pars$a_f2 * (1 - height / hmat_eff)))
}

## Read RA back out of the C++ model. There is no direct R binding for
## fraction_allocation_reproduction, but fecundity_dt = RA * dB/dt / (omega +
## a_f3) is exactly invertible, and both terms are exposed.
##
## This only works where net mass production is positive; below the whole-plant
## light compensation point compute_rates zeroes every rate and RA becomes
## unobservable. That is why L = 0 is not tested directly here -- a plant in
## the dark has no production to allocate. The L -> 0 limits (pow(0, 0) = 1
## giving plain FF16 when a_f5 = 0, and RA -> 0 when a_f5 > 0) follow from the
## closed form, which the parity test below pins down over the observable range.
ff16_ra_observed <- function(s, height, openness) {
  env <- Environment("FF16")
  env$set_fixed_environment(openness, 100)
  ind <- FF16_Individual(s)
  ind$set_state("height", height)
  ind$compute_rates(env)
  net_mass_production_dt <- ind$aux("net_mass_production_dt")
  expect_gt(net_mass_production_dt, 0)
  ind$rate("fecundity") * (s$pars$omega + s$pars$a_f3) / net_mass_production_dt
}

## Copy-back idiom: nested list access returns a copy, so modify it and assign
## it back. Assigning straight into `s$pars$a_f4` silently does nothing.
ff16_strategy_with <- function(a_f4 = 0, a_f5 = 0, a_f1 = 1, a_f2 = 50) {
  s <- FF16_Strategy()
  pars <- s$pars
  pars$a_f4 <- a_f4
  pars$a_f5 <- a_f5
  pars$a_f1 <- a_f1
  pars$a_f2 <- a_f2
  s$pars <- pars
  s
}

test_that("the copy-back idiom actually reaches the strategy", {
  s <- ff16_strategy_with(a_f4 = 0.5, a_f5 = 2, a_f1 = 0.6, a_f2 = 10)
  expect_equal(s$pars$a_f4, 0.5)
  expect_equal(s$pars$a_f5, 2)
  expect_equal(s$pars$a_f1, 0.6)
  expect_equal(s$pars$a_f2, 10)
})

test_that("a_f4 = a_f5 = 0 recovers FF16 bit-for-bit", {

  p0 <- scm_base_parameters("FF16")
  env <- Environment("FF16")
  ctrl <- Control()

  m_default <- trait_matrix(0.0825, "lma")
  ## Setting the new parameters explicitly also exercises them as traits --
  ## FF16_hyperpar passes through columns it does not generate, which is how
  ## the ESS sweeps will vary them.
  m_zero <- cbind(m_default, a_f4 = 0, a_f5 = 0)

  p_default <- add_strategies(p0, m_default, hyperpar = FF16_hyperpar,
                              birth_rate = list(20))
  p_zero <- add_strategies(p0, m_zero, hyperpar = FF16_hyperpar,
                           birth_rate = list(20))

  out_default <- run_scm(p_default, env, ctrl)
  out_zero <- run_scm(p_zero, env, ctrl)

  ## Bit-for-bit, not merely close.
  expect_identical(out_zero$offspring_production,
                   out_default$offspring_production)
  expect_identical(out_zero$ode_times, out_default$ode_times)

  ## ... and still on the pinned regression target.
  expect_equal(out_zero$offspring_production, 16.8895016, tolerance = 1e-4)
  expect_equal(out_zero$ode_times[c(10, 100)], c(0.000070, 4.215899),
               tolerance = 1e-5)
})

test_that("RA is independent of a_f4 and a_f5 in full sun", {
  height <- 17
  reference <- ff16_ra_observed(ff16_strategy_with(), height, 1.0)
  expect_equal(reference,
               ff16_ra_expected(FF16_Strategy()$pars, height, 1.0))

  for (a_f4 in c(0, 0.5, 2)) {
    for (a_f5 in c(0, 0.5, 2)) {
      s <- ff16_strategy_with(a_f4, a_f5)
      ## Exact: L = 1 makes L^a_f5 == 1 and hmat_eff == hmat for any values.
      expect_identical(ff16_ra_observed(s, height, 1.0), reference)
    }
  }
})

test_that("shading lowers RA when a_f4 or a_f5 is positive", {
  height <- 17
  ## Above the whole-plant light compensation point at this height, so
  ## production stays positive and RA remains readable.
  openness <- c(1.0, 0.9, 0.8, 0.7, 0.6)

  ra_of <- function(s) {
    vapply(openness, function(l) ff16_ra_observed(s, height, l), numeric(1))
  }

  ## Neither parameter set: RA depends on height alone.
  ra_flat <- ra_of(ff16_strategy_with())
  expect_equal(ra_flat, rep(ra_flat[[1]], length(openness)))

  ## Shade delay alone: RA falls monotonically as light falls.
  ra_delay <- ra_of(ff16_strategy_with(a_f4 = 0.5))
  expect_true(all(diff(ra_delay) < 0))
  expect_identical(ra_delay[[1]], ra_flat[[1]])

  ## Shade cap alone: likewise.
  ra_cap <- ra_of(ff16_strategy_with(a_f5 = 2))
  expect_true(all(diff(ra_cap) < 0))
  expect_identical(ra_cap[[1]], ra_flat[[1]])

  ## The shade cap acts as a pure multiplier L^a_f5 on the full-sun value.
  expect_equal(ra_cap / ra_flat, openness^2)
})

test_that("RA matches the closed-form reaction norm", {
  ## a_f1 and a_f2 are varied too: downstream code (regnans' fixed_RA) sets
  ## them, so the reaction norm has to compose with non-default values rather
  ## than only holding at FF16 defaults.
  grid <- expand.grid(height = c(14, 17, 20),
                      openness = c(1.0, 0.85, 0.7),
                      a_f4 = c(0, 0.5),
                      a_f5 = c(0, 1.5),
                      a_f1 = c(1, 0.6),
                      a_f2 = c(50, 10))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    s <- ff16_strategy_with(g$a_f4, g$a_f5, g$a_f1, g$a_f2)
    expect_equal(ff16_ra_observed(s, g$height, g$openness),
                 ff16_ra_expected(s$pars, g$height, g$openness),
                 tolerance = 1e-12,
                 info = paste(names(g), unlist(g), collapse = " "))
  }
})

test_that("shading still lowers RA when a_f1 and a_f2 are non-default", {
  height <- 17
  openness <- c(1.0, 0.9, 0.8, 0.7)
  ## a_f1 caps allocation below 1, a_f2 makes maturation gradual: the shape
  ## regnans' fixed_RA imposes.
  s <- ff16_strategy_with(a_f4 = 0.5, a_f5 = 1.5, a_f1 = 0.6, a_f2 = 10)
  ra <- vapply(openness, function(l) ff16_ra_observed(s, height, l), numeric(1))

  expect_true(all(diff(ra) < 0))
  ## Full sun is unaffected by the reaction norm, so it still sits at the
  ## height-only logistic with these a_f1/a_f2.
  expect_identical(ra[[1]],
                   ff16_ra_observed(ff16_strategy_with(a_f1 = 0.6, a_f2 = 10),
                                    height, 1.0))
  ## Allocation never exceeds the a_f1 ceiling.
  expect_true(all(ra <= 0.6))
})

test_that("Report generation", {

  p0 <- scm_base_parameters("FF16")
  env <- Environment("FF16")
  ctrl <- Control()
  
  p2 <- add_strategies(p0, trait_matrix(c(0.0825, 0.2625), "lma"), hyperpar = FF16_hyperpar, birth_rate = list(11.99177, 16.51006))

  # test report generation
  out <- run_scm(p2, env, ctrl, collect = TRUE)

  unlink("tmp", recursive = TRUE)
  expect_message(FF16_generate_stand_report(out, "tmp/tmp.html", overwrite = TRUE), "Report for FF16 stand saved at tmp/tmp.html")
  expect_true(file.exists("tmp/tmp.html"))

  # don't overwrite output if already exists 
  expect_message(FF16_generate_stand_report(out, "tmp/tmp.html", overwrite = FALSE), "Report for FF16 stand already exists at tmp/tmp.html")  
  expect_true(file.exists("tmp/tmp.html"))

  unlink("tmp", recursive = TRUE)

})
