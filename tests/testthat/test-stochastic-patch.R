
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

## ODE states contributed by the environment, and their values. Only TF24 carries
## any: five soil-water layers then five cumulative fluxes.
n_environment_ode_states <- c(FF16 = 0L, TF24 = 10L, K93 = 0L)

## Recompute the environment and every rate at the patch's current state, the way
## the solver does on each derivatives evaluation. Needed because the introduction
## path updates the competition profile without recomputing rates.
refresh_rates <- function(patch) {
  patch$set_ode_state(patch$ode_state, patch$ode_time)
  invisible(patch)
}

environment_ode_state <- function(patch, x) {
  if (n_environment_ode_states[[x]] == 0L) {
    return(numeric(0))
  }
  env <- patch$environment
  c(env$get_soil_water_state(), env$get_soil_water_state_cumulative_flux())
}

for (x in names(strategy_types)) {

  test_that("empty", {
  
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()))
    
    env <- Environment(x)
    ctrl <- Control()
    patch <- StochasticPatch(x, e)(p, env, ctrl)

    expect_inherits(patch, sprintf("StochasticPatch<%s,%s>",x,e))

    expect_equal(patch$size, 1)
    expect_equal(patch$height_max, 0.0)
    expect_equal(patch$get_area, 1.0)

    expect_equal(patch$compute_competition(0), 0.0)

    ## The ODE system is the individuals plus whatever state the environment
    ## carries, so an empty patch is the environment alone: nothing for FF16 and
    ## K93, TF24's five soil-water layers and four cumulative fluxes.
    expect_equal(patch$node_ode_size, 0)
    expect_equal(patch$ode_state, environment_ode_state(patch, x))
    expect_length(patch$ode_rates, n_environment_ode_states[[x]])
    expect_true(all(is.finite(patch$ode_rates)))

    sp <- patch$species
    expect_true(is.list(sp))
    expect_equal(length(sp), 1)
    expect_inherits(sp[[1]], sprintf("StochasticSpecies<%s,%s>",x,e))
    expect_equal(sp[[1]]$size, 0)
  })

  test_that("non empty", {
  
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()))
    
    env <- Environment(x)
    ctrl <- Control()
    patch <- StochasticPatch(x, e)(p, env, ctrl)
    cmp <- Individual(x, e)(p$strategies[[1]])
    
    expect_error(patch$introduce_new_node(0), "Invalid value")
    expect_error(patch$introduce_new_node(10), "out of bounds")

    # Deterministic insertion avoids flaky failures from establishment probability.
    patch$introduce_new_node_and_update(1)
    expect_gt(patch$height_max, 0.0)
    expect_equal(patch$height_max, cmp$state("height"))

    expect_equal(patch$deaths(), 0)

    ci <- patch$environment$light_availability$spline
    expect_equal(range(ci$x), c(0.0, cmp$state("height")))
    expect_equal(max(ci$y), 1.0)
    expect_lt(ci$y[[1]], 1.0)

    if (x == "FF16") {
      expect_true(all(patch$ode_rates > 0.0))
    }
  })


  test_that("change patch size", {
  
    env <- Environment(x)
    ctrl <- Control()
  
    e <- environment_types[[x]]
    p2 <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          patch_area= 2)
    patch2 <- StochasticPatch(x, e)(p2, env, ctrl)
    expect_equal(patch2$get_area, 2)
    expect_equal(p2$patch_area, 2)

    p10 <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          patch_area= 10)
    patch10 <- StochasticPatch(x, e)(p10, env, ctrl)
    expect_equal(p10$patch_area, 10)
    expect_equal(patch10$get_area, 10)
  })
}

## The stochastic patch's ODE system used to be the individuals alone, so an
## environment carrying state was held at its initial value for a whole run
## while the deterministic Patch integrated it. Only TF24 carries state, which
## is why this went unnoticed.
test_that("the environment's ODE state and rates are part of the system", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()))
    patch <- StochasticPatch(x, e)(p, Environment(x), Control())
    n_env <- n_environment_ode_states[[x]]

    expect_equal(patch$ode_size, n_env)

    patch$introduce_new_node_and_update(1)
    n_ind <- Individual(x, e)(p$strategies[[1]])$ode_size

    expect_equal(patch$node_ode_size, n_ind)
    expect_equal(patch$ode_size, n_ind + n_env)
    expect_length(patch$ode_state, n_ind + n_env)

    ## introduce_new_node_and_update() rebuilds the competition profile but does
    ## not recompute rates, so the environment's rates here are still the ones
    ## reset() left behind on an empty patch. Refresh through set_ode_state --
    ## the same door the solver uses -- so what follows tests computed rates.
    refresh_rates(patch)
    expect_length(patch$ode_rates, n_ind + n_env)
    expect_true(all(is.finite(patch$ode_rates)))

    ## The environment's states follow the individuals', as in Patch.
    expect_equal(tail(patch$ode_state, n_env), environment_ode_state(patch, x))
  }
})

## The two quantities StochasticPatch::compute_rates() introduces are the `/ area`
## normalisation of resource uptake and the summation over individuals -- each
## counted once, unlike Species, which integrates a density-weighted rate over the
## size distribution. Neither is pinned by asserting that uptake is positive and
## below infiltration: that holds for any positive uptake, so a factor error in
## either would land silently.
##
## They are testable together without reading an individual's uptake directly.
## Competition is leaf area per unit ground area, so N recruits on N m^2 present
## exactly the same competition profile as one on 1 m^2, and their physiology is
## therefore identical -- which the test asserts rather than assumes, by comparing
## the individuals' own rates. Holding the individuals fixed that way, per-area
## uptake must be invariant: the sum scales with N and the division by area
## divides it back out. Drop the division and the N = 4 case comes out 4x too
## high; count anything other than one uptake per living individual and it comes
## out wrong by that factor instead.
##
## The rate read is the cumulative root-uptake accumulator, whose rate is the
## per-area total across soil layers, with no positivity guard in the way. It is
## the second-to-last environment state: sum_pulse_runoff sits after it and is
## held at rate zero, so reading the last one would read a constant zero.
test_that("resource uptake is summed per individual and normalised by area", {
  patch_with <- function(patch_area, n_individuals) {
    p <- Parameters("TF24", "TF24_Env")(strategies = list(TF24_Strategy()),
                                       patch_area = patch_area)
    patch <- StochasticPatch("TF24", "TF24_Env")(p, Environment("TF24"), Control())
    for (i in seq_len(n_individuals)) {
      patch$introduce_new_node_and_update(1)
    }
    expect_equal(patch$species[[1]]$size, n_individuals)
    refresh_rates(patch)
    patch
  }
  uptake <- function(patch) tail(patch$ode_rates, 2)[[1]]
  first_individual_rates <- function(patch) {
    patch$species[[1]]$individual_at(1)$ode_rates
  }

  one <- patch_with(1, 1)
  ## Uptake has to be happening at all, or every comparison below is 0 == 0.
  expect_gt(uptake(one), 0)

  for (n in c(2, 4)) {
    dense <- patch_with(n, n)
    ## The premise: same competition per unit area, so the individuals are in the
    ## same environment and doing the same thing, to the last bit.
    expect_identical(first_individual_rates(dense), first_individual_rates(one))
    ## The conclusion: per-area uptake is unchanged, also to the last bit. Both
    ## the summation and the division have to be right for this to hold.
    expect_identical(uptake(dense), uptake(one))
  }

  ## And the sum really does reach every individual rather than stopping at the
  ## first: two on one square metre take close to twice as much. Not exactly
  ## twice, because TF24 averages light over the crown, so two individuals at the
  ## same height do shade each other a little even though neither is above the
  ## other -- which is why this direction is a bound rather than an equality.
  crowded <- uptake(patch_with(1, 2)) / uptake(one)
  expect_lt(crowded, 2)
  expect_gt(crowded, 1.99)
})

## An empty patch casts no shade, so compute_environment() discards the
## competition profile. It used to discard the environment's integrated state
## along with it, because TF24 restored its soil states from clear_environment().
## The solver evaluates derivatives on an empty patch whenever the schedule puts
## its first arrival some way into the run, so the state it integrated from was
## whatever the run started with.
test_that("an empty patch keeps the environment state the solver has set", {
  x <- "TF24"
  e <- environment_types[[x]]
  n_env <- n_environment_ode_states[[x]]
  patch <- StochasticPatch(x, e)(
    Parameters(x, e)(strategies=list(strategy_types[[x]]())),
    Environment(x), Control())

  expect_equal(patch$ode_size, n_env)
  moved <- environment_ode_state(patch, x)
  moved[[1]] <- moved[[1]] + 0.05

  patch$set_ode_state(moved, 1.0)
  expect_equal(environment_ode_state(patch, x), moved)
  expect_equal(patch$ode_state, moved)
})

## A recruit's strategy-specific initial states are seeded from its birth
## environment. The deterministic path does this in
## Node::compute_initial_conditions, before the first rates evaluation; the
## stochastic path used to skip it, so a TF24 seedling was born with an empty
## carbohydrate store. Mortality is not compared: establishment is a Bernoulli
## draw here rather than the -log(pr_estab) the deterministic node carries.
test_that("a stochastic recruit is born with its initial states seeded", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()))
    patch <- StochasticPatch(x, e)(p, Environment(x), Control())
    patch$introduce_new_node_and_update(1)

    cmp <- Individual(x, e)(p$strategies[[1]])
    cmp$set_initial_states(patch$environment)

    expect_equal(patch$species[[1]]$individual_at(1)$ode_state, cmp$ode_state)
  }
})

test_that("a TF24 recruit is born with reserves in its storage pool", {
  p <- Parameters("TF24", "TF24_Env")(strategies=list(TF24_Strategy()))
  patch <- StochasticPatch("TF24", "TF24_Env")(p, Environment("TF24"), Control())
  patch$introduce_new_node_and_update(1)

  ## a_st3 of storage capacity, which is what the deterministic node gets.
  node <- Node("TF24", "TF24_Env")(p$strategies[[1]])
  node$compute_initial_conditions(patch$environment, 1.0, 1.0)
  storage <- patch$species[[1]]$individual_at(1)$state("storage")

  expect_gt(storage, 0)
  expect_equal(storage, node$individual$state("storage"))
})
