context("StochasticPatch")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("empty", {
  for (x in names(strategy_types)) {
    context(sprintf("StochasticPatch-%s",x))

    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()))
    
    env <- make_environment(x)
    ctrl <- Control()
    patch <- StochasticPatch(x, e)(p, env, ctrl)

    expect_is(patch, sprintf("StochasticPatch<%s,%s>",x,e))

    expect_equal(patch$size, 1)
    expect_equal(patch$height_max, 0.0)
    expect_equal(patch$compute_competition(0), 0.0)
    expect_equal(patch$ode_state, numeric(0))
    expect_equal(patch$ode_rates, numeric(0))

    sp <- patch$species
    expect_true(is.list(sp))
    expect_equal(length(sp), 1)
    expect_is(sp[[1]], sprintf("StochasticSpecies<%s,%s>",x,e))
    expect_equal(sp[[1]]$size, 0)
  }
})

test_that("non empty", {
  for (x in names(strategy_types)) {
    context(sprintf("StochasticPatch-%s",x))

    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()))
    
    env <- make_environment(x)
    ctrl <- Control()
    patch <- StochasticPatch(x, e)(p, env, ctrl)
    cmp <- Individual(x, e)(p$strategies[[1]])
    
    patch

    expect_error(patch$introduce_new_node(0), "Invalid value")
    expect_error(patch$introduce_new_node(10), "out of bounds")

    expect_true(patch$introduce_new_node(1))
    expect_gt(patch$height_max, 0.0)
    expect_equal(patch$height_max, cmp$state("height"))

    expect_equal(patch$deaths(), 0)

    ci <- patch$environment$light_availability$spline
    expect_equal(range(ci$x), c(0.0, cmp$state("height")))
    expect_equal(max(ci$y), 1.0)
    expect_lt(ci$y[[1]], 1.0)

    if (x == "FF16") {
      expect_true(all(patch$ode_rates > 0.0))
    } else if (x == "FF16r") {
      expect_true(all(patch$ode_rates[-3] > 0.0))
      expect_equal(patch$ode_rates[[3]], 0)
    }
  }
})
