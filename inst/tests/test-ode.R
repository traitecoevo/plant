source("helper-tree.R")

context("ODE")

library(deSolve)
derivs.lorenz <- function(t, y, pars) {
  sigma <- pars[[1]]
  R <- pars[[2]]
  b <- pars[[3]]
  c(sigma * ( y[[2]] - y[[1]] ),
    R * y[[1]] - y[[2]] - y[[1]] * y[[3]],
    -b * y[[3]] + y[[1]] * y[[2]])
}
derivs.lorenz.d <- function(...)
  list(derivs.lorenz(...))

pars <- c(sigma=10.0,
          R=28.0,
          b=8.0 / 3.0)
t0 <- 0.0
tt <- seq(t0, t0+2, by=0.001)
y <- c(21, 21, 21)

lo <- new(Lorenz, pars[1], pars[2], pars[3])
expect_that(lo$derivs(t0, y),
            equals(derivs.lorenz(t0, y, pars)))

lo$set_state(y, t0)
expect_that(lo$state, is_identical_to(y))
expect_that(lo$time,  is_identical_to(t0))

## Try a fixed step, but one far too large:
lo$step_fixed(1)
expect_that(lo$state, is_identical_to(y))
expect_that(lo$time,  is_identical_to(t0))

dt <- 0.001
lo$step_fixed(dt)
## Force deSolve to take the same step.
cmp <- as.numeric(rk(y, c(0, dt), derivs.lorenz.d, pars,
                     method=rkMethod("rk45ck"), hini=dt, rtol=1,
                     atol=1)[2,-1])
expect_that(lo$state, equals(cmp, tolerance=1e-14))
expect_that(lo$time, is_identical_to(t0+dt))

## Variable step:
lo$set_state(y, t0)
lo$step()
expect_that(lo$time, is_greater_than(t0))
expect_that(identical(lo$state, y), is_false())

## Run:
lo$set_state(y, t0)
lo$advance(tt[2])
expect_that(lo$time, is_identical_to(tt[2]))
y.cmp <- lo$state

ans <- t(lo$run(tt, y))
expect_that(ans[1,], is_identical_to(y.cmp))

expect_that(nrow(ans), equals(length(tt)-1))

ans.d <- rk(y, tt, derivs.lorenz.d, pars,
            method=rkMethod("rk45ck"))[-1,-1,drop=FALSE]

expect_that(length(lo$times), is_greater_than(nrow(ans)))

expect_that(ans, equals(unname(ans.d), tolerance=1e-11))

y.curr <- lo$state
t.curr <- lo$time
t.next <- t.curr + 0.001
lo$step_to(t.next)
expect_that(lo$time, is_identical_to(t.next))
y.next <- lo$state
lo$set_state(y.curr, t.curr)
lo$advance(t.next)
expect_that(lo$state, equals(y.next))

## Rerun the analysis so that we can check the fixed spacing code:
ans <- t(lo$run(tt, y))
## At shis point, the steps vary so that they're slighty more tightly
## spaced at first, and relax down to being the same as the sample
## frequency.
t.steps <- lo$times

lo$set_state(y, t0)
lo$advance_fixed(t.steps)
expect_that(lo$time, is_identical_to(t.steps[[length(t.steps)]]))
expect_that(lo$state, is_identical_to(ans[nrow(ans),]))

## With the R Ode verison:

obj <- new(OdeR, derivs.lorenz, new.env(), pars)

expect_that(obj$derivs(t0, y),
            is_identical_to(derivs.lorenz(t0, y, pars)))

obj$set_state(y, t0)
expect_that(obj$state, is_identical_to(y))
expect_that(obj$time,  is_identical_to(t0))

## Try a fixed step, but one far too large:
obj$step_fixed(1)
expect_that(obj$state, is_identical_to(y))
expect_that(obj$time,  is_identical_to(t0))

dt <- 0.001
obj$step_fixed(dt)
## Force deSolve to take the same step.
cmp <- as.numeric(rk(y, c(0, dt), derivs.lorenz.d, pars,
                     method=rkMethod("rk45ck"), hini=dt, rtol=1,
                     atol=1)[2,-1])
expect_that(obj$state, equals(cmp, tolerance=1e-14))
expect_that(obj$time, is_identical_to(t0+dt))

## Variable step:
obj$set_state(y, t0)
obj$step()
expect_that(obj$time, is_greater_than(t0))
expect_that(identical(obj$state, y), is_false())

## Run:
obj$reset()
obj$set_state(y, t0)
obj$advance(tt[2])
expect_that(obj$time, is_identical_to(tt[2]))
y.cmp <- obj$state

ans <- t(obj$run(tt, y))
expect_that(ans[1,], is_identical_to(y.cmp))

expect_that(nrow(ans), equals(length(tt)-1))

ans.d <- rk(y, tt, derivs.lorenz.d, pars,
            method=rkMethod("rk45ck"))[-1,-1,drop=FALSE]

expect_that(ans, equals(unname(ans.d), tolerance=1e-11))

t.steps <- obj$times
obj$reset()
obj$set_state(y, t0)
obj$advance_fixed(t.steps)
expect_that(obj$time, is_identical_to(t.steps[[length(t.steps)]]))
expect_that(obj$state, is_identical_to(ans[nrow(ans),]))
expect_that(obj$times, is_identical_to(t.steps))

## Check that we are actually rewriting the time steps:
t.steps2 <- obj$times[1:10]
obj$set_state(y, t0)
obj$advance_fixed(t.steps2)
expect_that(obj$times, is_identical_to(t.steps2))

## And again, with control parameters:
ctrl <- new(Control)
ctrl.ode <- ctrl$ode_control
obj.ctrl <- new(OdeR, derivs.lorenz, new.env(), pars, ctrl.ode)

## This will be improved by not() in the next testthat.
differ <- !identical(obj$control$parameters,
                     obj.ctrl$control$parameters)
expect_that(differ, is_true())

## Run to make sure we've not broken anything obvious:
obj.ctrl$set_state(y, t0)
obj.ctrl$advance(tt[2])
expect_that(obj.ctrl$time, is_identical_to(tt[2]))
y.cmp <- obj.ctrl$state

ans <- t(obj.ctrl$run(tt, y))
expect_that(ans[1,], is_identical_to(y.cmp))

test_that("Infinite times cause errors, rather than infinite loops", {
  expect_that(obj$advance(Inf), throws_error())
  expect_that(obj$step_to(Inf), throws_error())
  expect_that(obj$step_fixed(Inf), throws_error())
})
