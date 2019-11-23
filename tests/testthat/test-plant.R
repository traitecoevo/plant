

strategy_types <- get_list_of_strategy_types()

for (x in names(strategy_types)) {

  context(sprintf("Plant-%s",x))

  test_that("Reference comparison", {
    s <- strategy_types[[x]]()
    pl <- Plant(x)(s)

    expect_is(pl, sprintf("Plant<%s>",x))
    # expect_is(pp, sprintf("PlantPlus<%s>",x))
    # expect_identical(pp$strategy, s)
    expect_identical(pl$strategy, s)

    ## Expected initial conditions
    h0 <- 10
    pl$set_state("height", h0)

    for(v in pl$ode_names[-1]) {
      expect_identical(pl$state(v), 0.0)
      expect_true(is.na(pl$rate(v)))
    }

    expect_identical(pl$state("height"), h0)

    ## Set and get functions behave identically
    pl$set_state("height", h0)

    expect_identical(pl$state("height"), h0)

    m0 <- 5
    pl$set_state("mortality", m0)

    expect_identical(pl$state("mortality"), m0)

    f0 <- 8
    pl$set_state("fecundity", f0)

    expect_identical(pl$state("fecundity"), f0)

    expect_error(pl$state("not a trait name"))

    ## Compare internals
    vars_pl <- pl$internals

    expect_is(vars_pl, "Internals")

    # variable_names <- c("area_leaf", "height", "mortality", "fecundity", "area_heartwood", "mass_heartwood")
    # rate_names <- paste0(setdiff(variable_names, "area_leaf"), "_dt")

    # expect_true(all(c(variable_names, rate_names) %in% names(vars_pl)))

    ## Compute the vital rates and compare them
    env <- test_environment(h0)
    light_env <- attr(env, "light_env") # underlying function

    pl$compute_rates(env)
    # pp$compute_rates(env)

    # vars_pp <- pp$internals
    vars_pl <- pl$internals

    # expect_equal(vars_pl[variable_names], vars_pp[variable_names])
    # expect_equal(vars_pl[rate_names], vars_pp[rate_names])

    ## Area_leaf_above
    # for (h in seq(0, h0, length.out=10)) {
    #   expect_identical(pl$compute_competition(h), pp$compute_competition(h))
    # }

    ## Germination_probability
    # expect_identical(pl$germination_probability(env), pp$germination_probability(env))

    ## ode_system
    # expect_identical(pl$ode_size, 5)
    
    # expect_identical(pl$ode_names, ode_names)

    ## ode_state
    expect_equal(pl$ode_state, c(h0, m0, f0, 0, 0))

    expect_equal(pl$strategy_name, x)
  })

  test_that("stochastic support", {

    s <- strategy_types[[x]]()
    p <- Plant(x)(s)

    expect_equal(p$state("mortality"), 0.0)
    expect_equal(p$mortality_probability, 0.0)

    p$set_state("mortality", pi)
    expect_equal(p$mortality_probability, 1 - exp(-pi))

    p$set_state("mortality", 10000000)
    expect_equal(p$mortality_probability, 1)

    p$reset_mortality()
    expect_equal(p$state("mortality"), 0.0)
    expect_equal(p$mortality_probability, 0.0)
  })

}
