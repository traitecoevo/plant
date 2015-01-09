if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("Strategy")

test_that("Defaults", {
  expected <- list(
    B1     = 0.306,
    B4     = 1.71,
    B5     = 0,
    B6     = 0,
    B7     = 1,
    Y      = 0.7,
    a1     = 5.44,
    a3     = 0.07,
    b      = 0.17,
    c_Rb   = 8024,
    c_Rl   = 21000,
    c_Rr   = 217,
    c_Rs   = 4012,
    c_acc  = 3.0*3.8e-5,
    c_bio  = 0.0245,
    c_d0   = 0.01,
    c_d1   = 0.0,
    c_d2   = 5.5,
    c_d3   = 20,
    c_p1   = 150.36,
    c_p2   = 0.19,
    c_r1   = 1,
    c_r2   = 50,
    c_s0   = 0.1,
    eta    = 12,
    hmat   = 16.5958691,
    hmat_0 = 16.5958691,
    k_b    = 0.2,
    k_l0   = 0.4565855,
    k_r    = 1,
    k_s0   = 0.2,
    lma    = 0.1978791,
    lma_0  = 0.1978791,
    n_area = 0.00187,
    n_area_0 = 1.87e-3,
    rho    = 608,
    rho_0  = 608,
    s      = 3.8e-5,
    s_0    = 3.8e-5,
    theta  = 4669,
    control = Control())

  keys <- sort(names(expected))

  s <- Strategy()
  expect_that(s, is_a("Strategy"))

  expect_that(sort(names(s)), is_identical_to(keys))
  expect_that(unclass(s)[keys], is_identical_to(expected[keys]))
})
