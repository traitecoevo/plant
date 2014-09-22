## Slow approach to equilibrium.
library(tree)
library(nleqslv)
library(rootSolve)
library(parallel)
source("equilibrium_fun.R")

## Here's a situation where our equilibrium finding totally falls apart.
time_disturbance <- 2.5
pars <- list(c_r1=0.5, c_r2=0, B5=0, c_d1=1, B6=1)

p <- make_pars(pars, time_disturbance)
add_strategy(p, list(rho=5.1))

res <- equilibrium_seed_rain(p)

approach <- t(sapply(attr(res, "progress"), "[[", "seed_rain"))
r <- range(approach)
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach, pch=19, cex=.5, type="o")

## # Case 1: single strategy, can exist at equilibrium:
f <- make_equilibrium_runner(p)

## This nails it:
tol <- p$control$parameters$equilibrium_eps
maxit <- p$control$parameters$equilibrium_nsteps
control <- list(xtol=tol, ftol=tol, maxit=maxit)

sol <- nleqslv(p$seed_rain, make_target(f), global="none", control=control)
out <- f(sol$x)
diff(as.numeric(out))

## Or with multiroot:
sol_m <- multiroot(make_target(f), p$seed_rain, positive=TRUE)

## # Case 2: single strategy, cannot exist at equilibrium:
p2 <- make_pars(pars, time_disturbance)
add_strategy(p2, list(rho=4))
f2 <- make_equilibrium_runner(p2)

## This one gets stuck, proposing a really negative value that we
## truncate and that's OK
sol2 <- nleqslv(p2$seed_rain, make_target(f2), global="none", control=control)

## multiroot actually does quite well here:
sol2_m <- multiroot(make_target(f2), p2$seed_rain, positive=TRUE)

## # Case 3: multiple strategies that displace each other:
p3 <- make_pars(pars, time_disturbance)
add_strategy(p3, list(rho=5.1))
add_strategy(p3, list(rho=5.102))
p3$seed_rain <- c(sol$x, 1)
f3 <- make_equilibrium_runner(p3)

## This gets into a bit of a tangle, suggesting values that don't
## really make much sense.  It makes really weird excursions because
## of the discontinuity at zero.
sol3 <- nleqslv(p3$seed_rain, make_target(f3), global="none", control=control)

## This gets tangled, then jumps to (0,0)
sol3_m <- multiroot(make_target(f3), p3$seed_rain, positive=TRUE)

## # Hail Mary idea:

## If we've introduced a mutant, we've already checked that it has
## positive fitness so we can exclude the zero solution for it.
allow_zero <- c(TRUE, FALSE)
g3 <- make_target2(f3, allow_zero)
g3(p3$seed_rain)

sol3 <- nleqslv(p3$seed_rain, g3, global="none", control=control)
sol3_m <- multiroot(g3, p3$seed_rain, positive=TRUE)
