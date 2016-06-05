context("PlantPlus")

strategy_types <- get_list_of_strategy_types()

test_that("Seed bits", {
  for (x in names(strategy_types)) {
    ## Seed stuff:
    s <- strategy_types[[x]]()
    seed <- PlantPlus(x)(s)
    env <- test_environment(10) # high enough
    light_env <- attr(env, "light_env") # underlying function

    ## Check that our root-finding succeeded and the leaf mass is correct:
    expect_equal(seed$internals[["mass_live"]], s$omega, tolerance=1e-7)

    cmp <- try(make_reference_plant(x), silent=TRUE)
    if (!inherits(cmp, "try-error")) {
      ## Check that the height at birth is correct.  These answers are
      ## actually quite different, which could come from the root finding?
      expect_equal(seed$height, cmp$height.at.birth(cmp$traits), tolerance=1e-4)

      ## Then, check the germination probabilities in the current light
      ## environment:
      expect_equal(seed$germination_probability(env),
                  cmp$germination.probability(cmp$traits, light_env),
                         tolerance=1e-5)
    }
  }
})

## TODO: Missing here: all the plant growing stuff.  Move that
## elsewhere.
test_that("Assimilation over distribution", {
  for (x in names(strategy_types)) {
    s1 <- strategy_types[[x]]()
    p1 <- PlantPlus(x)(s1)

    c2 <- Control(plant_assimilation_over_distribution=TRUE)
    s2 <- strategy_types[[x]](control=c2)
    p2 <- PlantPlus(x)(s2)

    p1$height <- 10.0
    p2$height <- p1$height
    env <- test_environment(p1$height)

    p1$compute_vars_phys(env)
    p2$compute_vars_phys(env)
    p1_vars <- p1$internals
    p2_vars <- p2$internals

    ## Result is similar but not identical:
    expect_equal(p2_vars, p1_vars, tolerance=1e-7)
    expect_not_identical(p2_vars, p1_vars)
  }
})

test_that("Non-adaptive assimilation integration works", {
  for (x in names(strategy_types)) {
    c1 <- Control(plant_assimilation_adaptive=TRUE,
                  plant_assimilation_over_distribution=TRUE)
    s1 <- strategy_types[[x]](control=c1)
    p1 <- PlantPlus(x)(s1)

    c2 <- Control(plant_assimilation_adaptive=FALSE,
                  plant_assimilation_over_distribution=TRUE)
    s2 <- strategy_types[[x]](control=c2)
    p2 <- PlantPlus(x)(s2)

    p1$height <- 10.0
    p2$height <- p1$height
    env <- test_environment(p1$height)

    p1$compute_vars_phys(env)
    p2$compute_vars_phys(env)
    p1_vars <- p1$internals
    p2_vars <- p2$internals

    ## Result is similar but not identical:
    expect_equal(p2_vars[["assimilation"]], p1_vars[["assimilation"]], tolerance=1e-3)
    expect_not_identical(p2_vars, p1_vars)
  }
})

test_that("Ode interface", {
  for (x in names(strategy_types)) {
    p <- PlantPlus(x)(strategy_types[[x]]())
    expect_equal(p$ode_size, 5)
    expect_equal(p$ode_state,
                 c(p$height, p$mortality, p$fecundity,
                         p$area_heartwood, p$mass_heartwood))

    env <- test_environment(p$height * 10)
    p$compute_vars_phys(env)
    p$compute_vars_growth() # NOTE: Compute immediately *after* vars_phys
    expect_equal(p$ode_state,
                 c(p$height, p$mortality, p$fecundity,
                         p$area_heartwood, p$mass_heartwood))
    vars <- as.list(p$internals)
    expect_equal(p$ode_rates,
                 c(vars$height_dt,
                         vars$mortality_dt,
                         vars$fecundity_dt,
                         vars$area_heartwood_dt,
                         vars$mass_heartwood_dt))

    state_new <- c(p$height * 2, runif(p$ode_size - 1L))
    p$ode_state <- state_new
    expect_identical(p$ode_state, state_new)

    expect_identical(p$ode_names,
                    c("height", "mortality", "fecundity",
                                  "area_heartwood", "mass_heartwood"))
  }
})

test_that("conversions", {
  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    p <- Plant(x)(s)
    p$height <- 10
    env <- test_environment(p$height)

    pp1 <- plant_to_plant_plus(p, NULL)
    expect_is(pp1, sprintf("PlantPlus<%s>",x))
    expect_equal(pp1$internals$assimilation, NA_real_)
    expect_equal(pp1$height, p$height)

    pp2 <- plant_to_plant_plus(p, env)
    expect_gt(pp2$internals$assimilation, 0)
    pp1$compute_vars_phys(env)
    pp1$compute_vars_growth()
    expect_equal(pp1$internals, pp2$internals)

    p2 <- pp1$to_plant()
    expect_is(p2, sprintf("Plant<%s>",x))
    expect_equal(p2$internals$height_dt, NA_real_)
    expect_equal(p2$internals, p$internals)
  }
})

test_that("regression", {
  s <- FF16_Strategy(
         hmat = 30.0,
         a_f1 = 0.8,
         a_f2 = 20,
         a_l1   = 2.17,
         a_l2   = 0.546,
         k_l  = 0.4565855 / 3,
         lma  = 0.06879341)
  pl <- FF16_PlantPlus(s)
  runner <- OdeRunner("PlantRunner")(PlantRunner(pl, fixed_environment(1)))
  runner$advance(5)
  d0 <- oderunner_plant_size(runner)[["diameter_stem"]]
  y0 <- runner$state
  t0 <- runner$time
  runner$step()
  runner$set_state(y0, t0)
  d1 <- oderunner_plant_size(runner)[["diameter_stem"]]
  expect_equal(d0, d1)
})
