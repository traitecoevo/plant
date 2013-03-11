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
expect_that(lo$get_state(), is_identical_to(y))
expect_that(lo$get_time(),  is_identical_to(t0))

## Try a fixed step, but one far too large:
lo$step_fixed(1)
expect_that(lo$get_state(), is_identical_to(y))
expect_that(lo$get_time(),  is_identical_to(t0))

dt <- 0.001
lo$step_fixed(dt)
## Force deSolve to take the same step.
cmp <- as.numeric(rk(y, c(0, dt), derivs.lorenz.d, pars,
                     method=rkMethod("rk45ck"), hini=dt, rtol=1,
                     atol=1)[2,-1])
expect_that(lo$get_state(), equals(cmp, tolerance=1e-14))
expect_that(lo$get_time(), is_identical_to(t0+dt))

## Variable step:
lo$set_state(y, t0)
lo$step()
expect_that(lo$get_time() > t0, is_true())
expect_that(identical(lo$get_state(), y), is_false())

## Run:
lo$set_state(y, t0)
lo$advance(tt[2])
expect_that(lo$get_time(), is_identical_to(tt[2]))
y.cmp <- lo$get_state()

ans <- t(lo$run(tt, y))
expect_that(ans[1,], is_identical_to(y.cmp))

expect_that(nrow(ans), equals(length(tt)-1))

ans.d <- rk(y, tt, derivs.lorenz.d, pars,
            method=rkMethod("rk45ck"))[-1,-1,drop=FALSE]

expect_that(ans, equals(unname(ans.d), tolerance=1e-11))

## With the R Ode verison:

obj <- new(ROde, derivs.lorenz, new.env(), pars)

expect_that(obj$derivs(t0, y),
            is_identical_to(derivs.lorenz(t0, y, pars)))

obj$set_state(y, t0)
expect_that(obj$get_state(), is_identical_to(y))
expect_that(obj$get_time(),  is_identical_to(t0))

## Try a fixed step, but one far too large:
obj$step_fixed(1)
expect_that(obj$get_state(), is_identical_to(y))
expect_that(obj$get_time(),  is_identical_to(t0))

dt <- 0.001
obj$step_fixed(dt)
## Force deSolve to take the same step.
cmp <- as.numeric(rk(y, c(0, dt), derivs.lorenz.d, pars,
                     method=rkMethod("rk45ck"), hini=dt, rtol=1,
                     atol=1)[2,-1])
expect_that(obj$get_state(), equals(cmp, tolerance=1e-14))
expect_that(obj$get_time(), is_identical_to(t0+dt))

## Variable step:
obj$set_state(y, t0)
obj$step()
expect_that(obj$get_time() > t0, is_true())
expect_that(identical(obj$get_state(), y), is_false())

## Run:
obj$set_state(y, t0)
obj$advance(tt[2])
expect_that(obj$get_time(), is_identical_to(tt[2]))
y.cmp <- obj$get_state()

ans <- t(obj$run(tt, y))
expect_that(ans[1,], is_identical_to(y.cmp))

expect_that(nrow(ans), equals(length(tt)-1))

ans.d <- rk(y, tt, derivs.lorenz.d, pars,
            method=rkMethod("rk45ck"))[-1,-1,drop=FALSE]

expect_that(ans, equals(unname(ans.d), tolerance=1e-11))
