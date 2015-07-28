context("Strategy-FFdev")

test_that("Defaults", {
  expected <- list(
    B1     = 0.306,
    Pi_0   = 0.25,
    Y      = 0.7,
    a1     = 5.44,
    a3     = 0.07,
    b      = 0.17,
    c_Rb   = 8024 / 608,
    c_Rl   = 39.27 / 0.1978791,
    c_Rr   = 217,
    c_Rs   = 4012/608,
    c_acc  = 3.0*3.8e-5,
    c_bio  = 0.0245,
    c_d0   = 0.01,
    c_d2   = 5.5,
    c_d3   = 20,
    c_p1   = 150.36,
    c_p2   = 0.19,
    c_r1   = 1,
    c_r3   = 2,
    c_s0   = 0.1,
    eta    = 12,
    hmat   = 16,
    k_b    = 0.2,
    k_l   = 0.4565855,
    k_r    = 1,
    k_s   = 0.2,
    lma    = 0.1978791,
    rho    = 608,
    mass_seed = 3.8e-5,
    theta  = 4669,
    control = Control())

  keys <- sort(names(expected))

  s <- FFdev_Strategy()
  expect_that(s, is_a("FFdev_Strategy"))

  expect_that(sort(names(s)), is_identical_to(keys))
  expect_that(unclass(s)[keys], is_identical_to(expected[keys]))
})

test_that("Reference comparison", {
  s <- FFdev_Strategy()
  p <- FFdev_PlantPlus(s)

  expect_that(p$strategy, is_identical_to(s))

  ## Set the height to something (here 10)
  h0 <- 10
  p$height <- h0

  vars <- p$internals

  expect_that(vars[["height"]],
              is_identical_to(h0))


  expect_that(p$height,    is_identical_to(vars[["height"]]))
  expect_that(p$area_leaf, is_identical_to(vars[["area_leaf"]]))

  # Reproductive allocation functions

  env <- test_environment(1)
  light_env <- attr(env, "light_env") # underlying function

  seed <- PlantPlus("FFdev")(s)

  h <- c(seed$height, s[["hmat"]],  s[["hmat"]] + s[["c_r3"]],  s[["hmat"]] + 3 * s[["c_r3"]])
  r_out <- c(0, 0, 0.5, 0.75)

  for (i in seq_along(h)) {
    p$height <- h[i]

    ## Compute the physiological variables and retrieve them.
    p$compute_vars_phys(env)
    p$compute_vars_growth()

    expect_that(p$internals[["fraction_allocation_reproduction"]],
      is_identical_to(r_out[i]))
  }

})
