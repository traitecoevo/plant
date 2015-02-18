context("Plant2")

test_that("Reference comparison", {
  s <- Strategy()
  p1 <- Plant(s)
  p2 <- Plant2(s)

  expect_that(p2$strategy, is_identical_to(s))
  expect_that(p2, is_a("Plant2"))

  h0 <- 10
  p2$height <- h0
  vars2 <- p2$internals
  expect_that(all(is.na(vars2[c("height_dt", "mortality_dt",
                                "fecundity_dt")])), is_true())

  for(v in c("mortality", "fecundity")) {
      expect_that(vars2[[v]], equals(0, tolerance=1e-7))
  }

  # Set and get functions behave identically
  p1$height <- h0
  p2$height <- h0

  vars1 <- p1$internals
  vars2 <- p2$internals

  expect_that(vars2, is_a("Plant2_internals"))

  variable.names <- c("area_leaf", "height", "height_dt", "mortality",
                      "mortality_dt","fecundity","fecundity_dt")

  expect_that(all(names(vars2) %in% variable.names), is_true())
  expect_that(all(names(vars2) %in% names(vars1)), is_true())

  for(v in variable.names) {
     expect_that(vars2[[v]], equals(vars1[[v]], tolerance=1e-7))
  }

  # Compute the vital rates and compare them

  env <- test_environment(h0)
  light_env <- attr(env, "light_env") # underlying function

  ## Compute the vital rates and compare them
  p1$compute_vars_phys(env)
  p2$compute_vars_phys(env)

  vars1 <- p1$internals
  vars2 <- p2$internals

  expect_that(vars1[["reproduction_dt"]],
              equals(vars2[["reproduction_dt"]], tolerance=1e-7))
  expect_that(vars1[["height_dt"]],
              equals(vars2[["height_dt"]], tolerance=1e-7))
  expect_that(vars1[["mortality_dt"]],
              equals(vars2[["mortality_dt"]], tolerance=1e-7))
})
