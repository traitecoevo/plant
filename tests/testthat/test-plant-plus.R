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
    expect_that(seed$internals[["mass_live"]],
                equals(s$mass_seed, tolerance=1e-7))

    cmp <- try(make_reference_plant(x), silent=TRUE)
    if (!inherits(cmp, "try-error")) {
      ## Check that the height at birth is correct.  These answers are
      ## actually quite different, which could come from the root finding?
      expect_that(seed$height,
                  equals(cmp$height.at.birth(cmp$traits), tolerance=1e-4))

      ## Then, check the germination probabilities in the current light
      ## environment:
      expect_that(seed$germination_probability(env),
                  equals(cmp$germination.probability(cmp$traits, light_env),
                         tolerance=1e-5))
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
    expect_that(p2_vars, equals(p1_vars, tolerance=1e-7))
    expect_that(p2_vars, not(is_identical_to(p1_vars)))
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
    expect_that(p2_vars[["assimilation"]],
                equals(p1_vars[["assimilation"]], tolerance=1e-3))
    expect_that(p2_vars, not(is_identical_to(p1_vars)))
  }
})

test_that("Ode interface", {
  for (x in names(strategy_types)) {
    p <- PlantPlus(x)(strategy_types[[x]]())
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
  }
})

test_that("conversions", {
  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    p <- Plant(x)(s)
    p$height <- 10
    env <- test_environment(p$height)

    pp1 <- plant_to_plant_plus(p, NULL)
    expect_that(pp1, is_a(sprintf("PlantPlus<%s>",x)))
    expect_that(pp1$internals$assimilation, equals(NA_real_))
    expect_that(pp1$height, equals(p$height))

    pp2 <- plant_to_plant_plus(p, env)
    expect_that(pp2$internals$assimilation, is_more_than(0))
    pp1$compute_vars_phys(env)
    pp1$compute_vars_growth()
    expect_that(pp1$internals, equals(pp2$internals))

    p2 <- pp1$to_plant()
    expect_that(p2, is_a(sprintf("Plant<%s>",x)))
    expect_that(p2$internals$height_dt, equals(NA_real_))
    expect_that(p2$internals, equals(p$internals))
  }
})

test_that("regression", {
  s <- FFW16_Strategy(
         hmat = 30.0,
         c_r1 = 0.8,
         c_r2 = 20,
         a1   = 2.17,
         B1   = 0.546,
         k_l  = 0.4565855 / 3,
         lma  = 0.06879341)
  pl <- FFW16_PlantPlus(s)
  runner <- OdeRunner("PlantRunner")(PlantRunner(pl, fixed_environment(1)))
  runner$advance(5)
  d0 <- oderunner_plant_size(runner)[["diameter_stem"]]
  y0 <- runner$state
  t0 <- runner$time
  runner$step()
  runner$set_state(y0, t0)
  d1 <- oderunner_plant_size(runner)[["diameter_stem"]]
  expect_that(d0, equals(d1))
})
