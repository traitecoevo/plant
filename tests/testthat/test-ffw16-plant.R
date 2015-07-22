context("FFW16_Plant")

test_that("Reference comparison", {
  s <- FFW16_Strategy()
  pl <- FFW16_Plant(s)
  pp <- FFW16_PlantPlus(s)

  expect_that(pl, is_a("Plant<FFW16>"))
  ## TODO: PlantPlus not yet templated (can it be?)
  expect_that(pp, is_a("FFW16_PlantPlus"))
  expect_that(pp$strategy, is_identical_to(s))
  expect_that(pl$strategy, is_identical_to(s))

  ## Expected initial conditions
  h0 <- 10
  pl$height <- h0
  vars_pl <- pl$internals
  expect_that(all(is.na(vars_pl[c("height_dt", "mortality_dt",
                                  "fecundity_dt",
                                  "area_heartwood_dt", "mass_heartwood_dt")])),
              is_true())
  for (v in c("mortality", "fecundity", "area_heartwood", "mass_heartwood")) {
    expect_that(vars_pl[[v]], is_identical_to(0.0))
  }
  expect_that(vars_pl$height, is_identical_to(h0))

  ## Set and get functions behave identically
  pl$height <- h0
  pp$height <- h0

  expect_that(pl$height, is_identical_to(h0))
  expect_that(pp$height, is_identical_to(h0))

  m0 <- 5
  pl$mortality <- m0
  pp$mortality <- m0

  expect_that(pl$mortality, is_identical_to(m0))
  expect_that(pp$mortality, is_identical_to(m0))

  f0 <- 8
  pl$fecundity <- f0
  pp$fecundity <- f0

  expect_that(pl$fecundity, is_identical_to(f0))
  expect_that(pp$fecundity, is_identical_to(f0))

  ## Compare internals
  vars_pl <- pl$internals
  vars_pp <- pp$internals

  expect_that(vars_pl, is_a("Plant_internals"))
  expect_that(vars_pp, is_a("FFW16_PlantPlus_internals"))

  variable_names <- c("area_leaf", "height", "mortality", "fecundity",
                      "area_heartwood", "mass_heartwood")
  rate_names <- paste0(setdiff(variable_names, "area_leaf"), "_dt")

  expect_that(all(c(variable_names, rate_names) %in% names(vars_pl)),
              is_true())
  expect_that(all(c(variable_names, rate_names) %in% names(vars_pp)),
              is_true())
  expect_that(all(names(vars_pl) %in% names(vars_pp)), is_true())

  expect_that(vars_pl[variable_names],
              is_identical_to(vars_pp[variable_names]))
  expect_that(vars_pl[rate_names],
              is_identical_to(vars_pp[rate_names]))

  ## Compute the vital rates and compare them
  env <- test_environment(h0)
  light_env <- attr(env, "light_env") # underlying function

  pl$compute_vars_phys(env)
  pp$compute_vars_phys(env)

  vars_pp <- pp$internals
  vars_pl <- pl$internals

  expect_that(vars_pl[variable_names],
              is_identical_to(vars_pp[variable_names]))
  expect_that(vars_pl[rate_names],
              is_identical_to(vars_pp[rate_names]))

  ## Area_leaf_above
  for (h in seq(0, h0, length.out=10)) {
    expect_that(pl$area_leaf_above(h),
                is_identical_to(pp$area_leaf_above(h)))
  }

  ## Germination_probability
  expect_that(pl$germination_probability(env),
              is_identical_to(pp$germination_probability(env)))

  ## ode_system
  expect_that(pl$ode_size, is_identical_to(5))
  expect_that(pp$ode_size, is_identical_to(5))

  ode_names <- c("height", "mortality", "fecundity",
                 "area_heartwood", "mass_heartwood")
  expect_that(pl$ode_names, is_identical_to(ode_names))
  expect_that(pp$ode_names, is_identical_to(ode_names))

  ## ode_rates
  expect_that(pl$ode_rates, equals(pp$ode_rates))

  ## ode_state
  expect_that(pl$ode_state, equals(c(h0, m0, f0, 0, 0)))
  expect_that(pl$ode_state, equals(pp$ode_state))
})

test_that("stochastic support", {
  s <- FFW16_Strategy()
  p <- FFW16_Plant(s)
  expect_that(p$mortality, equals(0.0))
  expect_that(p$mortality_probability, equals(0.0))

  p$mortality <- pi
  expect_that(p$mortality_probability, equals(1 - exp(-pi)))

  p$mortality <- 10000000
  expect_that(p$mortality_probability, equals(1))

  p$reset_mortality()
  expect_that(p$mortality, equals(0.0))
  expect_that(p$mortality_probability, equals(0.0))
})
