context("Strategy")

test_that("Defaults", {
  expected <- list(
    B1     = 0.306,
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
    c_d2   = 5.5,
    c_d3   = 20,
    c_p1   = 150.36,
    c_p2   = 0.19,
    c_r1   = 1,
    c_r2   = 50,
    c_s0   = 0.1,
    eta    = 12,
    hmat   = 16.5958691,
    k_b    = 0.2,
    k_l   = 0.4565855,
    k_r    = 1,
    k_s   = 0.2,
    lma    = 0.1978791,
    n_area = 0.00187,
    rho    = 608,
    s      = 3.8e-5,
    theta  = 4669,
    control = Control())

  keys <- sort(names(expected))

  s <- Strategy()
  expect_that(s, is_a("Strategy"))

  expect_that(sort(names(s)), is_identical_to(keys))
  expect_that(unclass(s)[keys], is_identical_to(expected[keys]))
})

test_that("Strategy parameters agree with reference model", {
  cmp <- make_reference_plant()
  cmp_pars <- cmp$get_parameters()

  s <- Strategy()

  ## Expect that all parameters in the R version are found in the C++
  ## version.
  expect_that(all(names(cmp_pars) %in% names(s)), is_true())

  ## And v.v., except for a few additions:
  extra <- "control"
  common <- setdiff(names(s), extra)
  expect_that(all(extra %in% names(s)), is_true())
  expect_that(all(common %in% names(cmp_pars)),
              is_true())

  ## The C++ version should have no NA values by this point.
  expect_that(any(sapply(s[common], is.na)),
              is_false())

  ## And neither should the R version.
  expect_that(any(sapply(cmp_pars, is.na)),
              is_false())

  ## And demand that all parameters agree.
  expect_that(s[names(cmp_pars)],
              is_identical_to(cmp_pars))
})
