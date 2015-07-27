context("FFdev_PlantPlus")

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

test_that("Seed bits", {
  ## Seed stuff:
  s <- FFdev_Strategy()
  seed <- FFdev_PlantPlus(s)
  env <- test_environment(10) # high enough
  light_env <- attr(env, "light_env") # underlying function

  ## Check that our root-finding succeeded and the leaf mass is correct:
  expect_that(seed$internals[["mass_live"]],
              equals(s$mass_seed, tolerance=1e-7))

})

## TODO: Missing here: all the plant growing stuff.  Move that
## elsewhere.

test_that("Assimilation over distribution", {
  s1 <- FFdev_Strategy()
  p1 <- FFdev_PlantPlus(s1)

  c2 <- Control(plant_assimilation_over_distribution=TRUE)
  s2 <- FFdev_Strategy(control=c2)
  p2 <- FFdev_PlantPlus(s2)

  p1$height <- 10.0
  p2$height <- p1$height
  env <- test_environment(p1$height)

  p1$compute_vars_phys(env)
  p2$compute_vars_phys(env)
  p1_vars <- p1$internals
  p2_vars <- p2$internals

  ## Result is similar but not identical:
  expect_that(p2_vars, equals(p1_vars, tolerance=1e-7))
  expect_that(p2_vars, not(is_identical_to(p1_vars)))
})

test_that("Non-adaptive assimilation integration works", {
  c1 <- Control(plant_assimilation_adaptive=TRUE,
                plant_assimilation_over_distribution=TRUE)
  s1 <- FFdev_Strategy(control=c1)
  p1 <- FFdev_PlantPlus(s1)

  c2 <- Control(plant_assimilation_adaptive=FALSE,
                plant_assimilation_over_distribution=TRUE)
  s2 <- FFdev_Strategy(control=c2)
  p2 <- FFdev_PlantPlus(s2)

  p1$height <- 10.0
  p2$height <- p1$height
  env <- test_environment(p1$height)

  p1$compute_vars_phys(env)
  p2$compute_vars_phys(env)
  p1_vars <- p1$internals
  p2_vars <- p2$internals

  ## Result is similar but not identical:
  expect_that(p2_vars[["assimilation"]],
              equals(p1_vars[["assimilation"]], tolerance=1e-3))
  expect_that(p2_vars, not(is_identical_to(p1_vars)))
})

test_that("Ode interface", {
  p <- FFdev_PlantPlus(FFdev_Strategy())
  expect_that(p$ode_size, equals(5))
  expect_that(p$ode_state,
              equals(c(p$height, p$mortality, p$fecundity,
                       p$area_heartwood, p$mass_heartwood)))

  env <- test_environment(p$height * 10)
  p$compute_vars_phys(env)
  p$compute_vars_growth() # NOTE: Compute immediately *after* vars_phys
  expect_that(p$ode_state,
              equals(c(p$height, p$mortality, p$fecundity,
                       p$area_heartwood, p$mass_heartwood)))
  vars <- as.list(p$internals)
  expect_that(p$ode_rates,
              equals(c(vars$height_dt,
                       vars$mortality_dt,
                       vars$fecundity_dt,
                       vars$area_heartwood_dt,
                       vars$mass_heartwood_dt)))

  state_new <- c(p$height * 2, runif(p$ode_size - 1L))
  p$ode_state <- state_new
  expect_that(p$ode_state, is_identical_to(state_new))

  expect_that(p$ode_names,
              is_identical_to(c("height", "mortality", "fecundity",
                                "area_heartwood", "mass_heartwood")))
})
