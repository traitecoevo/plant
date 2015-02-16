context("Lorenz (basic ODE)")

library(deSolve)
derivs_lorenz <- function(y, pars) {
  sigma <- pars[[1]]
  R <- pars[[2]]
  b <- pars[[3]]
  c(sigma * ( y[[2]] - y[[1]] ),
    R * y[[1]] - y[[2]] - y[[1]] * y[[3]],
    -b * y[[3]] + y[[1]] * y[[2]])
}

derivs_lorenz_d <- function(t, y, pars) {
  list(derivs_lorenz(y, pars))
}

test_that("Basic Lorenz object works", {
  pars <- c(sigma=10.0,
            R=28.0,
            b=8.0 / 3.0)
  t0 <- 0.0
  tt <- seq(t0, t0+2, by=0.001)
  y <- c(21, 21, 21)

  lo <- Lorenz(pars[[1]], pars[[2]], pars[[3]])
  ## Check that everything is sane to start off:
  expect_that(lo$ode_size, equals(3))
  expect_that(lo$pars, is_identical_to(unname(pars)))
  expect_that(lo$ode_state, is_identical_to(rep(0.0, 3)))
  expect_that(lo$ode_rates,
              is_identical_to(derivs_lorenz(lo$ode_state, pars)))
  expect_that(lo$ode_time, is_null())

  ## Then set the state:
  lo$ode_state <- y
  expect_that(lo$ode_state, is_identical_to(y))
  expect_that(lo$ode_rates, is_identical_to(derivs_lorenz(y, pars)))

  y2 <- runif(3)
  lo$ode_state <- y2
  expect_that(lo$ode_state, is_identical_to(y2))
  expect_that(lo$ode_rates, is_identical_to(derivs_lorenz(y2, pars)))
})

## Then, get the ode runner working.
test_that("Ode runner behaves", {
  pars <- c(sigma=10.0, R=28.0, b=8.0 / 3.0)
  y <- rep(21, 3)
  lo <- Lorenz(pars[[1]], pars[[2]], pars[[3]])
  lo$ode_state <- y
  derivs_lorenz(y, pars)

  sys <- OdeRunner("Lorenz")(lo)
  expect_that(sys, is_a("OdeRunner"))
  expect_that(sys, is_a("OdeRunner<Lorenz>"))

  expect_that(sys$time,  is_identical_to(0.0))
  expect_that(sys$state, is_identical_to(y))
  expect_that(sys$times, is_identical_to(0.0))
  expect_that(sys$object, is_a("Lorenz"))

  lo2 <- sys$object

  sys$step()
  expect_that(sys$time, is_more_than(0.0))
  expect_that(all(sys$state != y), is_true())
  expect_that(sys$object$ode_state, is_identical_to(sys$state))
  expect_that(lo2$ode_state, is_identical_to(y))

  t1 <- 1.0
  sys$advance(t1)

  ## State has changed:
  expect_that(all(sys$state != y), is_true())
  ## But not in these objects:
  expect_that(lo$ode_state,  is_identical_to(y))
  expect_that(lo2$ode_state, is_identical_to(y))

  times <- sys$times
  expect_that(first(times), is_identical_to(0.0))
  expect_that(last(times),  is_identical_to(t1))
  ## Object does not store time:
  expect_that(sys$object$ode_time, is_null())

  ## Run it again:
  sys2 <- OdeRunner("Lorenz")(lo)
  sys2$advance_fixed(times)
  expect_that(sys2$times, is_identical_to(times))
  expect_that(sys2$state, equals(sys$state, tolerance=1e-13))
})

test_that("OdeR interface", {
  pars <- c(sigma=10.0, R=28.0, b=8.0 / 3.0)
  y <- rep(21, 3)
  lo <- Lorenz(pars[[1]], pars[[2]], pars[[3]])
  lo$ode_state <- c(21, 21, 21)

  sys <- test_ode_make_system(lo)
  expect_that(sys, is_a("OdeR"))
  sol <- test_ode_make_solver(sys)
  expect_that(sol, is_a("OdeRunner"))
  expect_that(sol, is_a("OdeRunner<OdeR>"))

  expect_that(sol$state, is_identical_to(y))
  expect_that(sol$time,  is_identical_to(0.0))

  sol2 <- OdeRunner("Lorenz")(lo)

  sol$advance(pi)
  sol2$advance(pi)
  expect_that(sol$times, is_identical_to(sol2$times))
  expect_that(sol$state, is_identical_to(sol2$state))

  ## This *has* updated the original values.
  expect_that(lo$ode_state, is_identical_to(sol2$state))

  lo$ode_state <- y
  expect_that(lo$ode_state, not(is_identical_to(sol2$state)))
  expect_that(sol2$set_state_from_system(),
              throws_error("Time does not match previous"))
})
