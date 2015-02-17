## TODO: Test add_seeds(vector<double>)

context("Patch")

test_that("Basics", {
  ## TODO: This is something that needs validating: the seed_rain and
  ## is_resident vectors must be the right length.
  p <- Parameters(strategies=list(Strategy()),
                  seed_rain=pi/2,
                  is_resident=TRUE)
  patch <- Patch(p)
  cmp <- Cohort(p$strategies[[1]])

  expect_that(patch$size, equals(1))
  expect_that(patch$height_max, is_identical_to(cmp$height))
  expect_that(patch$parameters, equals(p))

  expect_that(patch$environment, is_a("Environment"))
  expect_that(patch$environment$time, is_identical_to(0.0))

  expect_that(length(patch$species), equals(1))
  expect_that(patch$species[[1]], is_a("Species"))

  expect_that(patch$ode_size, equals(0))
  expect_that(patch$ode_state, is_identical_to(numeric(0)))
  expect_that(patch$ode_rates, is_identical_to(numeric(0)))

  ## Empty light environment:
  patch$compute_light_environment()
  expect_that(patch$leaf_area_above(0), is_identical_to(0))

  expect_that(patch$add_seed(0), throws_error("Invalid value"))
  expect_that(patch$add_seed(2), throws_error("out of bounds"))

  patch$add_seed(1)
  expect_that(patch$ode_size, equals(4))

  ## Then pull this out:
  cmp$compute_initial_conditions(patch$environment)

  expect_that(patch$ode_state,
              is_identical_to(cmp$ode_state))
  expect_that(patch$ode_rates,
              is_identical_to(cmp$ode_rates))

  y <- patch$ode_state
  patch$set_ode_state(y, 0)
  expect_that(patch$ode_state, is_identical_to(y))

  ## NOTE: These should be identical, but are merely equal...
  expect_that(patch$derivs(y, 0),
              equals(cmp$ode_rates))

  ## solver <- solver_from_ode_target(patch, p$control$ode_control)
  ## solver$step()
  ## patch$add_seed(1)
  ## expect_that(patch$ode_size,
  ##             equals(cmp$ode_size * patch$n_individuals))

  patch$reset()
  expect_that(patch$ode_size, equals(0))
  expect_that(patch$environment$time, is_identical_to(0.0))

  t <- patch$environment$time # do via environment only?

  ## patch$add_seed(1)
  ## h <- patch$height[[1]]
  ## while (patch$time < 25) {
  ##   solver$step()
  ##   t <- c(t, patch$time)
  ##   h <- c(h, patch$height[[1]])
  ## }

  ## TODO: This is not really a test, but we need to look at this and
  ## see if it makes any sense at all.
  ## if (interactive()) {
  ##   plot(t, h, type="l")
  ##   plot(patch$environment$light_environment$xy, type="l")
  ## }

  ## patch$reset()
  ## patch$add_seed(1)
  ## solver <- solver_from_ode_target(patch, p$control$ode_control)

  ## tt <- seq(0, 25, length.out=26)
  ## hh <- patch$height[[1]]
  ## for (ti in tt[-1]) {
  ##   solver$advance(ti)
  ##   hh <- c(hh, patch$height[[1]])
  ## }

  ## if (interactive()) {
  ##   plot(t, h, type="l")
  ##   points(tt, hh)
  ## }

  ## expect_that(hh, equals(spline(t, h, xout=tt)$y, tolerance=1e-7))

  ## test_that("OK at end of sequence", {
  ##   expect_that(patch$time, is_identical_to(tt[[length(tt)]]))
  ##   solver$advance(tt[[length(tt)]])
  ##   expect_that(patch$time, is_identical_to(tt[[length(tt)]]))
  ##   expect_that(solver$advance(tt[[length(tt)]] - 1e-8),
  ##               throws_error())
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
  ##   patch$compute_vars_phys() # require because we just added seed
  ##   state <- patch$state

  ##   patch2 <- new(PatchCohortTop, patch$parameters)
  ##   expect_that(patch2$state <- state, throws_error())
  ##   patch2$force_state(state)
  ##   expect_that(patch2$state, is_identical_to(state))

  ##   ## Check some things that depend on state make sense:
  ##   expect_that(patch2$environment$light_environment$xy,
  ##               is_identical_to(patch$environment$light_environment$xy))
  ##   expect_that(patch2$time, is_identical_to(patch$time))
  ##   expect_that(patch2$ode_state, is_identical_to(patch$ode_state))
  ##   expect_that(patch2$height, is_identical_to(patch$height))
  ##   expect_that(patch2$ode_rates, is_identical_to(patch$ode_rates))
  ## })
})
