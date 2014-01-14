## Adaptive integration causes problems with the smoothness of the
## seed rain in/seed rain out function.
library(tree)
library(parallel)

## Try to establish what the equilibrium seed rain is.  First, a
## non-adaptive set of parameters
p.n <- new(Parameters)
p.n$add_strategy(new(Strategy))
p.n$seed_rain <- 1.1                       # Starting rain.
p.n$set_parameters(list(patch_area=1.0))   # See issue #13
p.n$set_control_parameters(fast.control()) # A bit faster

## And a set that is identical except for the adaptive assimilation is
## turned on:
p.a <- p.n$copy()
p.a$set_control_parameters(list(plant_assimilation_adaptive=TRUE))

## And a seed set of cohort introduction times for schedule building.
t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)

## Functions to run the model with an input seed rain, parameter set
## and times.  The second one will take the times as the seed point to
## build a new cohort schedule.
run <- function(seed_rain.in, p, times) {
  p$seed_rain <- seed_rain.in
  run.ebt(p, schedule.from.times(times))$fitness(1)
}

run.new.schedule <- function(w, p, times, build.args=list()) {
  p$seed_rain <- w
  build.args <- modifyList(list(nsteps=20, eps=1e-3, verbose=FALSE),
                           build.args)
  res <- build.schedule(p, times, build.args$nsteps, build.args$eps,
                        progress=FALSE, verbose=build.args$verbose)
  attr(res, "seed_rain")[,"out"]
}

## # 1: Approach to equilibrium:
res.n <- equilibrium.seed.rain(p.n, times0, 10, progress=TRUE,
                               build.args=list(verbose=TRUE))
res.a <- equilibrium.seed.rain(p.a, times0, 10, progress=TRUE,
                               build.args=list(verbose=TRUE))

approach.n <- t(sapply(attr(res.n, "progress"), "[[", "seed_rain"))
approach.a <- t(sapply(attr(res.a, "progress"), "[[", "seed_rain"))

cols <- c(n="black", a="red")
w.hat <- unname(res.n[["seed_rain"]][,"out"])

## From a distance, these both hone in nicely on the equilibrium, and
## rapidly, too.
r <- range(approach.n, approach.a)
plot(approach.n, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach.n, col=cols[["n"]], pch=19, cex=.5, type="o")
cobweb(approach.a, col=cols[["a"]], pch=19, cex=.5, type="o")

## Here is the painful lack of convergence for the
## with-adaptive-integration case (red); it's converging to a *region*
## but moving essentially stochastically within that region.
r <- w.hat + c(-1, 1) * 0.5
plot(approach.n, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach.n, col=cols[["n"]], pch=19, cex=.5, type="o")
cobweb(approach.a, col=cols[["a"]], pch=19, cex=.5, type="o")

## # 2: Near the equilibrium point:

## Then, in the vinicity of the root we should look at what the curve
## actually looks like, without adaptive refinement.
dw <- 2 # range of input to vary (plus and minus this many seeds)
seed_rain.in  <- seq(w.hat - dw, w.hat + dw, length=31)

times1 <- res.n[["times"]]
seed_rain.out.n <- unlist(mclapply(seed_rain.in, run, p.n, times1))
seed_rain.out.a <- unlist(mclapply(seed_rain.in, run, p.a, times1))

fit.n <- lm(seed_rain.out.n ~ seed_rain.in)
fit.a <- lm(seed_rain.out.a ~ seed_rain.in)

matplot(seed_rain.in, cbind(seed_rain.out.n, seed_rain.out.a),
        xlab="Incoming seed rain", ylab="Outgoing seed rain",
        las=1, pch=1, col=cols)
abline(0, 1, lty=2, col="grey")
abline(fit.n, lty=2, col=cols[["n"]])
abline(fit.a, lty=2, col=cols[["a"]])
cobweb(approach.n, col=cols[["n"]])
cobweb(approach.a, col=cols[["a"]])

## # 3: Global function shape
seed_rain.in.global <- seq(1, w.hat + 10, length=51)

seed_rain.out.global.n <-
  unlist(mclapply(seed_rain.in.global, run.new.schedule, p.n, times0))
seed_rain.out.global.a <-
  unlist(mclapply(seed_rain.in.global, run.new.schedule, p.a, times0))

matplot(seed_rain.in.global,
        cbind(seed_rain.out.global.n, seed_rain.out.global.a),
        las=1, type="l", col=cols, lty=1,
        xlab="Incoming seed rain", ylab="Outgoing seed rain")
abline(0, 1, lty=2, col="grey")
abline(fit.n, lty=2, col=cols[["n"]])
abline(fit.a, lty=2, col=cols[["a"]])
cobweb(approach.n, col=cols[["n"]])
cobweb(approach.a, col=cols[["a"]])
