# Built from  tests/testthat/test-strategy-ff16.R on Fri Oct 30 11:43:30 2020 using the scaffolder, from the strategy:  FF16
context("Strategy-FF16bg")

test_that("Defaults", {
  expected <- list(
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
    k_2 = 0.0,
    recruitment_decay = 0,
    control = Control(),
    collect_all_auxiliary = FALSE,
    birth_rate_x = numeric(0), # empty
    birth_rate_y = c(1.0),
    is_variable_birth_rate = FALSE)

  keys <- sort(names(expected))

  s <- FF16bg_Strategy()
  expect_is(s, "FF16bg_Strategy")

  expect_identical(sort(names(s)), keys)
  expect_identical(unclass(s)[keys], expected[keys])
})

test_that("FF16bg collect_all_auxiliary option", {

  s <- FF16bg_Strategy()
  p <- FF16bg_Individual(s)
  expect_equal(p$aux_size, 2)
  expect_equal(length(p$internals$auxs), 2)
  expect_equal(p$aux_names, c(
    "competition_effect",
    "net_mass_production_dt"
  ))

  s <- FF16bg_Strategy(collect_all_auxiliary=TRUE)
  expect_true(s$collect_all_auxiliary)
  p <- FF16bg_Individual(s)
  expect_equal(p$aux_size, 3)
  expect_equal(length(p$internals$auxs), 3)
  expect_equal(p$aux_names, c(
    "competition_effect",
    "net_mass_production_dt",
    "area_sapwood"
  ))
})

test_that("Reference comparison", {
  s <- FF16bg_Strategy()
  p <- FF16bg_Individual(s)

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
  s <- FF16bg_Strategy()
  my_names <- FF16bg_Individual(s)$ode_names
  expect_identical(my_names[1:3], c("height", "mortality", "fecundity"))
})
test_that("FF16bg_Strategy hyper-parameterisation", {
  s <- FF16bg_Strategy()

  # lma
  lma <- c(0.1,1)
  ret <- FF16bg_hyperpar(trait_matrix(lma, "lma"), s)

  expect_true(all(c("lma", "k_l", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "lma"], lma)
  expect_equal(ret[, "k_l"], c(1.46678,0.028600), tolerance=1e-5)
  expect_equal(ret[, "r_l"], c(392.70, 39.27), tolerance=1e-5)

  ## This happens on Linux (and therefore on travis) due to numerical
  ## differences in the integration.
  if ("a_p1" %in% colnames(ret)) {
    a_p1 <- ret[, "a_p1"]
    expect_equal(length(unique(a_p1)), 1L)
    expect_equal(a_p1[[1]], s$a_p1, tolerance=1e-7)
  }

  # wood density
  rho <- c(200,300)
  ret <- FF16bg_hyperpar(trait_matrix(rho, "rho"), s)
  expect_true(all(c("rho", "r_s", "r_b") %in% colnames(ret)))
  expect_equal(ret[, "rho"], rho)
  expect_equal(ret[, "r_s"], c(20.06000,13.37333), tolerance=1e-5)
  expect_equal(ret[, "r_b"], 2*ret[, "r_s"])

  ## This happens on Linux (and therefore on travis) due to numerical
  ## differences in the integration.
  if ("a_p1" %in% colnames(ret)) {
    a_p1 <- ret[, "a_p1"]
    expect_equal(length(unique(a_p1)), 1L)
    expect_equal(a_p1[[1]], s$a_p1, tolerance=1e-7)
  }

  # narea
  narea <- c(0, 2E-3,2.3E-3)
  ret <- FF16bg_hyperpar(trait_matrix(narea, "narea"), s)
  expect_true(all(c("narea", "a_p1", "a_p2", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "narea"], narea)
  expect_equal(ret[, "r_l"], c(0, 212.2508, 244.0884), tolerance=1e-5)
  expect_equal(ret[, "a_p1"], c(0, 162.2592, 188.1549), tolerance=1e-5)
  expect_equal(ret[, "a_p2"], c(0, 0.220904, 0.259173), tolerance=1e-5)

  # seed mass
  omega <- 3.8e-5*c(1,2,3)
  ret <- FF16bg_hyperpar(trait_matrix(omega, "omega"), s)
  expect_true(all(c("omega", "a_f3") %in% colnames(ret)))
  expect_equal(ret[, "omega"], omega)
  expect_equal(ret[, "a_f3"], 3*omega)

  ## This happens on Linux (and therefore on travis) due to numerical
  ## differences in the integration.
  if ("a_p1" %in% colnames(ret)) {
    a_p1 <- ret[, "a_p1"]
    expect_equal(length(unique(a_p1)), 1L)
    expect_equal(a_p1[[1]], s$a_p1, tolerance=1e-7)
  }


  ## Empty trait matrix:
  ret <- FF16bg_hyperpar(trait_matrix(numeric(0), "lma"), s)
  expect_equal(ret, trait_matrix(numeric(0), "lma"))
})

test_that("narea calculation", {
  x <- c(1.38, 3.07, 2.94)
  p0 <- FF16bg_Parameters()
  m <- trait_matrix(x, "hmat")
  expect_silent(sl <- strategy_list(m, p0, FF16bg_hyperpar, birth_rate_list=1.0))

  cmp <- lapply(x, function(xi) strategy(trait_matrix(xi, "hmat"), p0, FF16bg_hyperpar, birth_rate_list=1.0))
  expect_equal(sl, cmp)
})

# integration test - runs a full patch metapopultaion
# the offspring production produced integrates all demographic behaviours

test_that("offspring production", {

  p0 <- scm_base_parameters("FF16")
  env <- make_environment("FF16")
  ctrl <- scm_base_control()

  # Should default to FF16
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16bg_hyperpar,
                          mutant = FALSE, birth_rate_list = list(20))

  out <- run_scm(p1, env, ctrl)
  expect_equal(out$offspring_production, 16.88946, tolerance=1e-5)
  expect_equal(out$ode_times[c(10, 100)], c(0.000070, 4.216055), tolerance=1e-5)

  # Increasing below-ground importance
  p2 <- expand_parameters(trait_matrix(0.05, "k_2"), p0, FF16bg_hyperpar,
                          mutant = FALSE, birth_rate_list = list(20))

  out <- run_scm(p2, env, ctrl)
  expect_equal(out$offspring_production, 48.5288, tolerance=1e-5)
  expect_equal(out$ode_times[c(10, 100)], c(0.000070, 4.768224), tolerance=1e-5)
})
