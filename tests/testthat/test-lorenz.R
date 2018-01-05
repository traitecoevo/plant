context("Lorenz (basic ODE)")

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
  expect_equal(lo$ode_size, 3)
  expect_identical(lo$pars, unname(pars))
  expect_identical(lo$ode_state, rep(0.0, 3))
  expect_identical(lo$ode_rates, derivs_lorenz(lo$ode_state, pars))
  expect_null(lo$ode_time)

  ## Then set the state:
  lo$ode_state <- y
  expect_identical(lo$ode_state, y)
  expect_identical(lo$ode_rates, derivs_lorenz(y, pars))

  y2 <- runif(3)
  lo$ode_state <- y2
  expect_identical(lo$ode_state, y2)
  expect_identical(lo$ode_rates, derivs_lorenz(y2, pars))
})

## Then, get the ode runner working.
test_that("Ode runner behaves", {
  pars <- c(sigma=10.0, R=28.0, b=8.0 / 3.0)
  y <- rep(21, 3)
  lo <- Lorenz(pars[[1]], pars[[2]], pars[[3]])
  lo$ode_state <- y
  derivs_lorenz(y, pars)

  sys <- OdeRunner("Lorenz")(lo)
  expect_is(sys, "OdeRunner")
  expect_is(sys, "OdeRunner<Lorenz>")

  expect_identical(sys$time, 0.0)
  expect_identical(sys$state, y)
  expect_identical(sys$times, 0.0)
  expect_is(sys$object, "Lorenz")

  lo2 <- sys$object

  sys$step()
  expect_gt(sys$time, 0.0)
  expect_true(all(sys$state != y))
  expect_identical(sys$object$ode_state, sys$state)
  expect_identical(lo2$ode_state, y)

  t1 <- 1.0
  sys$advance(t1)

  ## State has changed:
  expect_true(all(sys$state != y))
  ## But not in these objects:
  expect_identical(lo$ode_state, y)
  expect_identical(lo2$ode_state, y)

  times <- sys$times
  expect_identical(first(times), 0.0)
  expect_identical(last(times), t1)
  ## Object does not store time:
  expect_null(sys$object$ode_time)

  ## Run it again:
  sys2 <- OdeRunner("Lorenz")(lo)
  sys2$advance_fixed(times)
  expect_identical(sys2$times, times)
  expect_equal(sys2$state, sys$state, tolerance=1e-13)
})

test_that("OdeR interface", {
  pars <- c(sigma=10.0, R=28.0, b=8.0 / 3.0)
  y <- rep(21, 3)
  lo <- Lorenz(pars[[1]], pars[[2]], pars[[3]])
  lo$ode_state <- c(21, 21, 21)

  sys <- test_ode_make_system(lo)
  expect_is(sys, "OdeR")
  sol <- test_ode_make_solver(sys)
  expect_is(sol, "OdeRunner")
  expect_is(sol, "OdeRunner<OdeR>")

  expect_identical(sol$state, y)
  expect_identical(sol$time, 0.0)

  sol2 <- OdeRunner("Lorenz")(lo)

  sol$advance(pi)
  sol2$advance(pi)
  expect_identical(sol$times, sol2$times)
  expect_identical(sol$state, sol2$state)

  ## This *has* updated the original values.
  expect_identical(lo$ode_state, sol2$state)

  lo$ode_state <- y
  expect_false(identical(lo$ode_state, sol2$state))
  expect_error(sol2$set_state_from_system(), "Time does not match previous")
})
