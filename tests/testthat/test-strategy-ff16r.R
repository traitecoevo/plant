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

    control = Control())

  keys <- sort(names(expected))

  s <- FF16r_Strategy()
  expect_is(s, "FF16r_Strategy")

  expect_identical(sort(names(s)), keys)
  expect_identical(unclass(s)[keys], expected[keys])
})

test_that("Reference comparison", {
  s <- FF16r_Strategy()
  p <- FF16r_PlantPlus(s)

  expect_identical(p$strategy, s)

  ## Set the height to something (here 10)
  h0 <- 10
  p$height <- h0

  vars <- p$internals

  expect_identical(vars[["height"]], h0)


  expect_identical(p$height, vars[["height"]])
  expect_identical(p$area_leaf, vars[["area_leaf"]])

  # Reproductive allocation functions

  env <- test_environment(1)
  light_env <- attr(env, "light_env") # underlying function

  seed <- PlantPlus("FF16r")(s)

  h <- c(seed$height, s[["hmat"]],  s[["hmat"]] + s[["a_f2"]],  s[["hmat"]] + 3 * s[["a_f2"]])
  r_out <- c(0, 0, 0.5, 0.75)

  for (i in seq_along(h)) {
    p$height <- h[i]

    ## Compute the physiological variables and retrieve them.
    p$compute_vars_phys(env)
    p$compute_vars_growth()

    expect_identical(p$internals[["fraction_allocation_reproduction"]], r_out[i])
  }

})
