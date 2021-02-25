## TODO: The tests here really warrant splitting into different chunks
## - this was ported over from tree1 where the tests were loose in the
## file.

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  e <- environment_types[[x]]

  context(sprintf("Species-%s",x))

  test_that("Basics", {
    env <- test_environment(x, 3, seed_rain=1.0)
    s <- strategy_types[[x]]()
    sp <- Species(x, e)(s)
    seed <- Cohort(x, e)(s)
    plant <- Individual(x, e)(s)
    h0 <- seed$height

    expect_equal(sp$size, 0)
    expect_identical(sp$height_max, h0)
    expect_identical(sp$cohorts, list())
    expect_identical(sp$height, NULL)
    expect_identical(sp$log_densities, numeric(0))
    expect_identical(sp$competition_effects, numeric(0))
    expect_identical(sp$competition_effects_error(1.0), numeric(0))
    expect_equal(sp$ode_size, 0)
    expect_identical(sp$ode_state, numeric(0))
    expect_identical(sp$ode_rates, numeric(0))

    ## Causes initial conditions to be estimated:
    sp$compute_rates(env)
    seed$compute_initial_conditions(env)

    ## Internal and test seed report same values:
    expect_identical(sp$seed$rates, seed$rates)
    expect_identical(sp$seed$ode_state, seed$ode_state)

    sp$add_seed()
    expect_equal(sp$size, 1)

    cohorts <- sp$cohorts
    expect_is(cohorts, "list")
    expect_equal(length(cohorts), 1)
    expect_identical(cohorts[[1]]$rates, seed$rates)
    expect_equal(sp$heights, seed$height)
    expect_equal(sp$log_densities, seed$log_density)
    expect_equal(sp$competition_effects, seed$competition_effect)
    ## NOTE: Didn't check ode values

    ## Internal and test seed report same values:
    expect_identical(sp$seed$rates, seed$rates)

    expect_is(sp$cohort_at(1), sprintf("Cohort<%s,%s>",x,e))
    expect_identical(sp$cohort_at(1)$rates, cohorts[[1]]$rates)

    ## Not sure about this -- do we need more immediate access?
    expect_identical(sp$seed$plant$establishment_probability(env), plant$establishment_probability(env))

    expect_equal(sp$compute_competition(0), 0)

    sp$heights <- 1

    h <- 0
    xx <- c(sp$seed$height, sp$heights)
    y <- c(sp$seed$compute_competition(h),
           sp$cohort_at(1)$compute_competition(h))

    expect_identical(sp$compute_competition(h), trapezium(xx, y))

    ## Better tests: I want cases where:
    ## 1. empty: throws error
    sp$clear()

    ## Re-set up the initial conditions
    sp$compute_rates(env)

    expect_error(sp$cohort_at(1), "Index 1 out of bounds")
    expect_error(sp$cohort_at(0), "Invalid value for index")
  })

  ## 1: empty species (no cohorts) has no leaf area above any height:
  test_that("Empty species has no leaf area", {
    sp <- Species(x, e)(strategy_types[[x]]())
    expect_equal(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(10), 0)
    expect_equal(sp$compute_competition(Inf), 0)
  })

  ## 2: Cohort up against boundary has no leaf area:
  test_that("species with only boundary cohort no leaf area", {
    env <- test_environment(x, 3, seed_rain=1.0)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$add_seed()
    sp$compute_rates(env)
    expect_equal(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(10), 0)
    expect_equal(sp$compute_competition(Inf), 0)
  })

  cmp_compute_competition <- function(h, sp) {
    x <- c(sp$heights, sp$seed$height)
    y <- c(sapply(sp$cohorts, function(p) p$compute_competition(h)),
           sp$seed$compute_competition(h))
    trapezium(rev(x), rev(y))
  }

  ## 3: Single cohort; one round of trapezium:
  test_that("Leaf area sensible with one cohort", {
    env <- test_environment(x, 3, seed_rain=1.0)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$compute_rates(env)
    sp$add_seed()
    h_top <- sp$height_max * 4
    sp$heights <- h_top

    ## At base and top
    expect_gt(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(0), cmp_compute_competition(0, sp))

    expect_identical(sp$compute_competition(h_top), 0.0)

    ## Part way up (and above bottom seed boundary condition)
    expect_equal(sp$compute_competition(h_top * .5), cmp_compute_competition(h_top * .5, sp))

    ode_size <- Cohort(x, e)(strategy_types[[x]]())$ode_size
    ode_state <- sp$ode_state
    p <- sp$cohort_at(1)
    expect_equal(sp$ode_size, ode_size)
    expect_equal(length(ode_state), ode_size)
    expect_identical(ode_state, p$ode_state)
  })

  test_that("Leaf area sensible with two cohorts", {
    env <- test_environment(x, 3, seed_rain=1.0)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$compute_rates(env)
    sp$add_seed()
    h_top <- sp$height_max * 4
    sp$add_seed()
    sp$heights <- h_top * c(1, .6)

    ## At base and top
    expect_gt(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(0), cmp_compute_competition(0, sp))
    expect_equal(sp$compute_competition(h_top), 0)
    ## Part way up (below bottom cohort, above boundarty condition)
    expect_equal(sp$compute_competition(h_top * .5), cmp_compute_competition(h_top * .5, sp))
    ## Within the top pair (excluding the seed)
    expect_equal(sp$compute_competition(h_top * .8), cmp_compute_competition(h_top * .8, sp))

    ode_size <- Cohort(x, e)(strategy_types[[x]]())$ode_size
    ode_state <- sp$ode_state
    cohorts <- sp$cohorts
    expect_equal(sp$ode_size, ode_size * sp$size)
    expect_equal(length(ode_state), ode_size * sp$size)
    expect_identical(ode_state, unlist(lapply(cohorts, function(p) p$ode_state)))
  })

  test_that("Leaf area sensible with three cohorts", {
    env <- test_environment(x, 3, seed_rain=1.0)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$compute_rates(env)
    sp$add_seed()
    h_top <- sp$height_max * 4
    sp$add_seed()
    sp$add_seed()
    sp$heights <- h_top * c(1, .75, .6)

    ## At base and top
    expect_gt(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(0), cmp_compute_competition(0, sp))
    expect_equal(sp$compute_competition(h_top), 0)
    ## Part way up (below bottom cohort, above boundarty condition)
    expect_equal(sp$compute_competition(h_top * .5), cmp_compute_competition(h_top * .5, sp))
    ## Within the top pair (excluding the seed)
    expect_equal(sp$compute_competition(h_top * .8), cmp_compute_competition(h_top * .8, sp))

    cmp_competition_effect <- sapply(seq_len(sp$size),
                            function(i) sp$cohort_at(i)$competition_effect)
    expect_identical(sp$competition_effects, cmp_competition_effect)

    cmp    <- local_error_integration(sp$heights, cmp_competition_effect, 1.0)
    cmp_pi <- local_error_integration(sp$heights, cmp_competition_effect, pi)

    expect_identical(sp$competition_effects_error(), cmp)
    expect_identical(sp$competition_effects_error(1.0), cmp)
    expect_identical(sp$competition_effects_error(pi), cmp_pi)

    ode_size <- Cohort(x, e)(strategy_types[[x]]())$ode_size
    ode_state <- sp$ode_state
    cohorts <- sp$cohorts
    expect_equal(length(ode_state), ode_size * sp$size)
    expect_identical(ode_state, unlist(lapply(cohorts, function(p) p$ode_state)))
  })
}
