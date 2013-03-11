source("helper-tree.R")

context("ODE")

library(deSolve)
derivs.lorenz <- function(t, y, pars) {
  sigma <- pars[1]
  R <- pars[2]
  b <- pars[3]
  list(c(sigma * ( y[2] - y[1] ),
         R * y[1] - y[2] - y[1] * y[3],
         -b * y[3] + y[1] * y[2]))
}
pars <- c(sigma=10.0,
          R=28.0,
          b=8.0 / 3.0)
t0 <- 0.0
tt <- seq(t0, t0+2, by=0.001)
y <- c(21, 21, 21)
## ans <- rk(y, tt, derivs.lorenz, pars, method = rkMethod("rk45ck"))

## xx <- ans[,2]
## yy <- ans[,3]
## zz <- ans[,4]
## plt <- persp(matrix(1, 2, 2), border=NA, col=NA, xlim=range(xx),
##              ylim=range(yy), zlim=range(zz), theta=30,
##              phi=60, box=FALSE, xlab="x", ylab="y", zlab="z")
## lines(trans3d(xx, yy, zz, plt))

e <- new(Evolve)

e$set_state(y, t0)
expect_that(e$get_state(), is_identical_to(y))
expect_that(e$get_time(),  is_identical_to(t0))

## Check the derivatives function:
expect_that(e$derivs(),
            equals(unname(derivs.lorenz(t0, y, pars)[[1]])))

## Try a fixed step, but one far too large:
e$step_fixed(1)
expect_that(e$get_state(), is_identical_to(y))
expect_that(e$get_time(),  is_identical_to(t0))

dt <- 0.001
e$step_fixed(dt)
## Force deSolve to take the same step.
cmp <- as.numeric(rk(y, c(0, dt), derivs.lorenz, pars,
                     method=rkMethod("rk45ck"), hini=dt, rtol=1,
                     atol=1)[2,-1])
expect_that(e$get_state(), equals(cmp, tolerance=1e-14))
expect_that(e$get_time(), is_identical_to(t0+dt))

## Variable step:
e$set_state(y, t0)
e$step()
expect_that(e$get_time() > t0, is_true())
expect_that(identical(e$get_state(), y), is_false())

## Run:
e$set_state(y, t0)
e$advance(tt[2])
expect_that(e$get_time(), is_identical_to(tt[2]))
y.cmp <- e$get_state()

ans <- t(e$run(tt, y))
expect_that(ans[1,], is_identical_to(y.cmp))

expect_that(nrow(ans), equals(length(tt)-1))

ans.d <- rk(y, tt, derivs.lorenz, pars,
            method=rkMethod("rk45ck"))[-1,-1,drop=FALSE]

expect_that(ans, equals(unname(ans.d), tolerance=1e-11))
