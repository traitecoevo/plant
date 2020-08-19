# Built from  tests/testthat/test-strategy-ff16.R on Fri Jul  3 08:14:35 2020 using the scaffolder, from the strategy:  FF16
context("Strategy-FF16r")

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
    a_f2   = 2,
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
    control = Control(),
    collect_all_auxillary = FALSE)

  keys <- sort(names(expected))

  s <- FF16r_Strategy()
  expect_is(s, "FF16r_Strategy")

  expect_identical(sort(names(s)), keys)
  expect_identical(unclass(s)[keys], expected[keys])
})

test_that("FF16r collect_all_auxillary option", {

  s <- FF16r_Strategy()
  p <- FF16r_Individual(s)
  expect_equal(p$aux_size, 2)
  expect_equal(length(p$internals$auxs), 2)
  expect_equal(p$aux_names, c(
    "competition_effect",
    "net_mass_production_dt"
  ))

  s <- FF16r_Strategy(collect_all_auxillary=TRUE)
  expect_true(s$collect_all_auxillary)
  p <- FF16r_Individual(s)
  expect_equal(p$aux_size, 3)
  expect_equal(length(p$internals$auxs), 3)
  expect_equal(p$aux_names, c(
    "competition_effect",
    "net_mass_production_dt",
    "area_sapwood"
  ))
})

test_that("Reference comparison", {
  s <- FF16r_Strategy()
  p <- FF16r_Individual(s)

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
  s <- FF16r_Strategy()
  my_names <- FF16r_Individual(s)$ode_names
  expect_identical(my_names[1:3], c("height", "mortality", "fecundity"))
})
test_that("FF16r_Strategy hyper-parameterisation", {
  s <- FF16r_Strategy()

  # lma
  lma <- c(0.1,1)
  ret <- FF16r_hyperpar(trait_matrix(lma, "lma"), s)

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
  ret <- FF16r_hyperpar(trait_matrix(rho, "rho"), s)
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
  ret <- FF16r_hyperpar(trait_matrix(narea, "narea"), s)
  expect_true(all(c("narea", "a_p1", "a_p2", "r_l") %in% colnames(ret)))
  expect_equal(ret[, "narea"], narea)
  expect_equal(ret[, "r_l"], c(0, 212.2508, 244.0884), tolerance=1e-5)
  expect_equal(ret[, "a_p1"], c(0, 162.2592, 188.1549), tolerance=1e-5)
  expect_equal(ret[, "a_p2"], c(0, 0.220904, 0.259173), tolerance=1e-5)

  # seed mass
  omega <- 3.8e-5*c(1,2,3)
  ret <- FF16r_hyperpar(trait_matrix(omega, "omega"), s)
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
  ret <- FF16r_hyperpar(trait_matrix(numeric(0), "lma"), s)
  expect_equal(ret, trait_matrix(numeric(0), "lma"))
})

