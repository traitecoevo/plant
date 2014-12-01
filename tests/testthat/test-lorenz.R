if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

## What I need to provide here is a way to set up a list of
## non-interacting Lorenz attractors with something controlling them.
## http://headmyshoulder.github.io/odeint-v2/doc/boost_numeric_odeint/tutorial/self_expanding_lattices.html
##
## Once everything is set up, might be worth trying an adams bashfoth
## approach.  That doesn't work perfectly with resizing, but could
## save a lot of effort?
##
## Probably need something to check that the state of the system has
## actually changed.
##
## Exploit the FSAL
## http://headmyshoulder.github.io/odeint-v2/doc/boost_numeric_odeint/odeint_in_detail/steppers.html
##
## Also exploit the fact that we know dydt on the way in: there's a
## method for that!

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
  expect_equal(lo$size, 3)
  expect_identical(lo$pars, pars)
  expect_identical(lo$ode_values, rep(0.0, 3))
  expect_identical(lo$ode_rates, derivs_lorenz(lo$ode_values, pars))

  ## Then set the state:
  lo$ode_values <- y
  expect_identical(lo$ode_values, y)
  expect_identical(lo$ode_rates, derivs_lorenz(y, pars))
})

## Then, get the ode runner working.
test_that("Ode runner behaves", {
  pars <- c(sigma=10.0, R=28.0, b=8.0 / 3.0)
  y <- c(21, 21, 21)
  lo <- Lorenz(pars[[1]], pars[[2]], pars[[3]])
  lo$ode_values <- y
  
  sys <- OdeSystem("Lorenz")(lo)
  lo2 <- sys$obj
  expect_identical(lo2$ode_values, lo$ode_values)
  y2 <- runif(3)
  lo2$ode_values <- y2
  expect_identical(lo2$ode_values, y2)
  ## Unchanged:
  expect_identical(lo$ode_values, y)
  expect_identical(sys$obj$ode_values, y)

  ## OK, what I've not done yet is shown that we can make this step
  ## *correctly*, deal with time stepping issues, or use a given
  ## schedule.  All these things need doing.
  sys$do_step(0.001)
  sys$obj$ode_values
  sys$y
  sys$try_step(1)

  sys$advance(1, 0.01)
  sys$y

  res <- sys$advance_save(10, 0.001)
  ## pairs(t(res), panel=lines)
})
