context("Plant2")

test_that("Reference comparison", {
  s <- Strategy()
  p1 <- Plant(s)
  p2 <- Plant2(s)

  expect_that(p2$strategy, is_identical_to(s))
  expect_that(p2, is_a("Plant2"))

  h0 <- 10
  p1$height <- h0
  p2$height <- h0

  vars1 <- p1$internals
  vars2 <- p2$internals

  expect_that(vars2, is_a("Plant2_internals"))

  expect_that(p1$height, is_identical_to(p2$height))
  expect_that(vars2[["height"]], is_identical_to(vars1[["height"]]))
  expect_that(vars2[["leaf_area"]], is_identical_to(vars1[["leaf_area"]]))

  env <- test_environment(h0)
  light_env <- attr(env, "light_env") # underlying function

  ## Compute the vital rates and compare them
  p1$compute_vars_phys(env)
  p2$compute_vars_phys(env)

  vars1 <- p1$internals
  vars2 <- p2$internals

  expect_that(vars1[["fecundity_rate"]],
              equals(vars2[["fecundity_rate"]], tolerance=1e-7))
  expect_that(vars1[["height_growth_rate"]],
              equals(vars2[["height_rate"]], tolerance=1e-7))
  expect_that(vars1[["mortality_rate"]],
              equals(vars2[["mortality_rate"]], tolerance=1e-7))
})
