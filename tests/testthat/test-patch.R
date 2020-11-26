## TODO: Test add_seeds(vector<double>)

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  context(sprintf("Patch-%s",x))

  test_that(sprintf("Basics %s", x), {
    ## TODO: This is something that needs validating: the seed_rain and
    ## is_resident vectors must be the right length.

    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    plant <- Individual(x, e)(s)
    cohort <- Cohort(x, e)(s)

    p <- Parameters(x, e)(strategies=list(s),
                          seed_rain=pi/2,
                          is_resident=TRUE)

    patch <- Patch(x, e)(p)
    cmp <- Cohort(x, e)(p$strategies[[1]])

    expect_equal(patch$size, 1)
    expect_identical(patch$height_max, cmp$height)
    expect_equal(patch$parameters, p)

    # This doesn't hold now
    # expect_is(patch$environment, paste0(x, "_Environment"))
    expect_identical(patch$environment$time, 0.0)

    expect_equal(length(patch$species), 1)
    expect_is(patch$species[[1]], sprintf("Species<%s,%s>",x,e))

    expect_equal(patch$ode_size, 0)
    expect_identical(patch$ode_state, numeric(0))
    expect_identical(patch$ode_rates, numeric(0))

    ## Empty light environment:
    patch$compute_environment()
    expect_identical(patch$compute_competition(0), 0)

    expect_error(patch$add_seed(0), "Invalid value")
    expect_error(patch$add_seed(2), "out of bounds")

    ode_size <- Cohort(x, e)(s)$ode_size
    patch$add_seed(1)
    expect_equal(patch$ode_size, ode_size)

    ## Then pull this out:
    cmp$compute_initial_conditions(patch$environment)

    expect_identical(patch$ode_state, cmp$ode_state)
    expect_identical(patch$ode_rates, cmp$ode_rates)

    y <- patch$ode_state
    patch$set_ode_state(y, 0)
    expect_identical(patch$ode_state, y)

    ## NOTE: These should be identical, but are merely equal...
    expect_equal(patch$derivs(y, 0), cmp$ode_rates)

    ## solver <- solver_from_ode_target(patch, p$control$ode_control)
    ## solver$step()
    ## patch$add_seed(1)
    ## expect_equal(patch$ode_size,
    ##             cmp$ode_size * patch$n_individuals)

    patch$reset()
    expect_equal(patch$ode_size, 0)
    expect_identical(patch$environment$time, 0.0)

    t <- patch$environment$time # do via environment only?

    ## patch$add_seed(1)
    ## h <- patch$state("height")[[1]]
    ## while (patch$time < 25) {
    ##   solver$step()
    ##   t <- c(t, patch$time)
    ##   h <- c(h, patch$state("height")[[1]])
    ## }

    ## TODO: This is not really a test, but we need to look at this and
    ## see if it makes any sense at all.
    ## if (interactive()) {
    ##   plot(t, h, type="l")
    ##   plot(patch$environment$environment_interpolator$xy, type="l")
    ## }

    ## patch$reset()
    ## patch$add_seed(1)
    ## solver <- solver_from_ode_target(patch, p$control$ode_control)

    ## tt <- seq(0, 25, length.out=26)
    ## hh <- patch$state("height")[[1]]
    ## for (ti in tt[-1]) {
    ##   solver$advance(ti)
    ##   hh <- c(hh, patch$state("height")[[1]])
    ## }

    ## if (interactive()) {
    ##   plot(t, h, type="l")
    ##   points(tt, hh)
    ## }

    ## expect_equal(hh, spline(t, h, xout=tt)$y, tolerance=1e-7)

    ## test_that("OK at end of sequence", {
    ##   expect_identical(patch$time, tt[[length(tt)]])
    ##   solver$advance(tt[[length(tt)]])
    ##   expect_identical(patch$time, tt[[length(tt)]])
    ##   expect_error(solver$advance(tt[[length(tt)]] - 1e-8))
    ## })

    ## test_that("State get/set works", {
    ##   patch$reset()
    ##   patch$add_seed(1)
    ##   ode.control <- p$control$ode_control
    ##   ode.control$set_parameters(list(step_size_min = 1e-4))
    ##   solver <- solver_from_ode_target(patch, ode.control)
    ##   while (patch$time < 5) {
    ##     solver$step()
    ##     if (patch$time > patch$n_individuals) {
    ##       patch$add_seed(1)
    ##       solver <- solver_from_ode_target(patch, ode.control)
    ##     }
    ##   }
    ##   patch$compute_rates() # require because we just added seed
    ##   state <- patch$state

    ##   patch2 <- new(PatchCohortTop, patch$parameters)
    ##   expect_error(patch2$state <- state)
    ##   patch2$force_state(state)
    ##   expect_identical(patch2$state, state)

    ##   ## Check some things that depend on state make sense:
    ##   expect_identical(patch2$environment$environment_interpolator$xy)
    ##   expect_identical(patch2$time, patch$time)
    ##   expect_identical(patch2$ode_state, patch$ode_state)
    ##   expect_identical(patch2$height, patch$state("height"))
    ##   expect_identical(patch2$ode_rates, patch$ode_rates)
    ## })
  })
}
