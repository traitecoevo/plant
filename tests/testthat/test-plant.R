

strategy_types <- get_list_of_strategy_types()

for (x in names(strategy_types)) {

  context(sprintf("Plant-%s",x))

  test_that("Reference comparison", {
    s <- strategy_types[[x]]()
    pl <- Plant(x)(s)
    pp <- PlantPlus(x)(s)

    expect_is(pl, sprintf("Plant<%s>",x))
    expect_is(pp, sprintf("PlantPlus<%s>",x))
    expect_identical(pp$strategy, s)
    expect_identical(pl$strategy, s)

    ## Expected initial conditions
    h0 <- 10
    pl$height <- h0
    vars_pl <- pl$internals
    expect_true(all(is.na(vars_pl[c("height_dt", "mortality_dt",
                                    "fecundity_dt",
                                    "area_heartwood_dt", "mass_heartwood_dt")])))
    for (v in c("mortality", "fecundity", "area_heartwood", "mass_heartwood")) {
      expect_identical(vars_pl[[v]], 0.0)
    }
    expect_identical(vars_pl$height, h0)

    ## Set and get functions behave identically
    pl$height <- h0
    pp$height <- h0

    expect_identical(pl$height, h0)
    expect_identical(pp$height, h0)

    m0 <- 5
    pl$mortality <- m0
    pp$mortality <- m0

    expect_identical(pl$mortality, m0)
    expect_identical(pp$mortality, m0)

    f0 <- 8
    pl$fecundity <- f0
    pp$fecundity <- f0

    expect_identical(pl$fecundity, f0)
    expect_identical(pp$fecundity, f0)

    ## Compare internals
    vars_pl <- pl$internals
    vars_pp <- pp$internals

    expect_is(vars_pl, "Plant_internals")
    expect_is(vars_pp, "PlantPlus_internals")

    variable_names <- c("area_leaf", "height", "mortality", "fecundity",
                        "area_heartwood", "mass_heartwood")
    rate_names <- paste0(setdiff(variable_names, "area_leaf"), "_dt")

    expect_true(all(c(variable_names, rate_names) %in% names(vars_pl)))
    expect_true(all(c(variable_names, rate_names) %in% names(vars_pp)))
    expect_true(all(names(vars_pl) %in% names(vars_pp)))

    expect_identical(vars_pl[variable_names], vars_pp[variable_names])
    expect_identical(vars_pl[rate_names], vars_pp[rate_names])

    ## Compute the vital rates and compare them
    env <- test_environment(h0)
    light_env <- attr(env, "light_env") # underlying function

    pl$compute_vars_phys(env)
    pp$compute_vars_phys(env)

    vars_pp <- pp$internals
    vars_pl <- pl$internals

    expect_equal(vars_pl[variable_names], vars_pp[variable_names])
    expect_equal(vars_pl[rate_names], vars_pp[rate_names])

    ## Area_leaf_above
    for (h in seq(0, h0, length.out=10)) {
      expect_identical(pl$area_leaf_above(h), pp$area_leaf_above(h))
    }

    ## Germination_probability
    expect_identical(pl$germination_probability(env), pp$germination_probability(env))

    ## ode_system
    expect_identical(pl$ode_size, 5)
    expect_identical(pp$ode_size, 5)

    ode_names <- c("height", "mortality", "fecundity",
                   "area_heartwood", "mass_heartwood")
    expect_identical(pl$ode_names, ode_names)
    expect_identical(pp$ode_names, ode_names)

    ## ode_rates
    expect_equal(pl$ode_rates, pp$ode_rates)

    ## ode_state
    expect_equal(pl$ode_state, c(h0, m0, f0, 0, 0))
    expect_equal(pl$ode_state, pp$ode_state)
  })

  test_that("stochastic support", {

    s <- strategy_types[[x]]()
    p <- Plant(x)(s)

    expect_equal(p$mortality, 0.0)
    expect_equal(p$mortality_probability, 0.0)

    p$mortality <- pi
    expect_equal(p$mortality_probability, 1 - exp(-pi))

    p$mortality <- 10000000
    expect_equal(p$mortality_probability, 1)

    p$reset_mortality()
    expect_equal(p$mortality, 0.0)
    expect_equal(p$mortality_probability, 0.0)
  })

}
