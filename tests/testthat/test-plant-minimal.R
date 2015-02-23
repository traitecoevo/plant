context("PlantMinimal")

test_that("Reference comparison", {
  s <- Strategy()
  p1 <- Plant(s)
  p2 <- PlantMinimal(s)

  expect_that(p2, is_a("PlantMinimal"))
  expect_that(p2$strategy, is_identical_to(s))
  expect_that(p2$strategy, is_identical_to(p1$strategy))

  # Expected initial conditions
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

  expect_that(p2$height, is_identical_to(h0))
  expect_that(p1$height, is_identical_to(p2$height))

  m0 <- 5
  p1$mortality <- m0
  p2$mortality <- m0

  expect_that(p2$mortality, is_identical_to(m0))
  expect_that(p1$mortality, is_identical_to(p2$mortality))

  f0 <- 8
  p1$fecundity <- f0
  p2$fecundity <- f0

  expect_that(p2$fecundity, is_identical_to(f0))
  expect_that(p1$fecundity, is_identical_to(p2$fecundity))


  # Compare internals
  vars1 <- p1$internals
  vars2 <- p2$internals

  expect_that(vars2, is_a("PlantMinimal_internals"))

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

  p1$compute_vars_phys(env)
  p2$compute_vars_phys(env)

  vars1 <- p1$internals
  vars2 <- p2$internals

  for(v in variable.names) {
     expect_that(vars2[[v]], equals(vars1[[v]], tolerance=1e-7))
  }

  # Area_leaf_above
  for(h in seq(0, h0, length.out=10)){
    expect_that(p2$area_leaf_above(h), equals(p1$area_leaf_above(h),
                                              tolerance=1e-7))
  }

  # Germination_probability
  expect_that(p2$germination_probability(env),
              equals(p1$germination_probability(env), tolerance=1e-7))

  # ode_system
  expect_that(p2$ode_size, is_identical_to(3))
  expect_that(p2$ode_size, is_identical_to(p1$ode_size -2))

  ode_names <- c("height", "mortality", "fecundity")
  expect_that(p2$ode_names, is_identical_to(ode_names))
  expect_that(all(p2$ode_names %in% p1$ode_names), is_true())
  ii <- match(ode_names, p1$ode_names)

  # ode_rates
  expect_that(p2$ode_rates, equals(p1$ode_rates[ii]))

  # ode_state
  expect_that(p2$ode_state, equals(c(h0, m0, f0)))
  expect_that(p2$ode_state, equals(p1$ode_state[ii]))
})
