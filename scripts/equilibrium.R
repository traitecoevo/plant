library(tree)
library(parallel)

## Try to establish what the equilubrium seed rain is.
p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- 1.1                       # Starting rain.
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

times0 <- cohort.introduction.times(104)
schedule0 <- schedule.from.times(times0)

run <- function(seed_rain.in, p, schedule) {
  p$seed_rain <- seed_rain.in
  run.ebt(p, schedule)$fitnesses
}

run.new.schedule <- function(w, p, schedule, build.args=list()) {
  p$seed_rain <- w
  build.args <- modifyList(list(nsteps=20, eps=1e-3, verbose=FALSE),
                           build.args)
  res <- build.schedule(p, schedule, build.args$nsteps, build.args$eps,
                        progress=FALSE, verbose=build.args$verbose)
  attr(res, "seed_rain")[,"out"]
}

## # 1: Approach to equilibrium:
res <- equilibrium.seed.rain(p, schedule0, 10, progress=TRUE,
                             build.args=list(verbose=TRUE))
w.hat <- unname(res[["seed_rain"]][,"out"])

## Plot the approach:
approach <- t(sapply(attr(res, "progress"), "[[", "seed_rain"))

## From a distance, these both hone in nicely on the equilibrium, and
## rapidly, too.
r <- range(approach)
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach, pch=19, cex=.5, type="o")

## Zoom in on the last few points:
r <- w.hat + c(-1, 1) * 0.03
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach, pch=19, cex=.5, type="o")

## # 2: Near the equilibrium point:

## Then, in the vinicity of the root we should look at what the curve
## actually looks like, without adaptive refinement.
dw <- 2 # range of input to vary (plus and minus this many seeds)
seed_rain.in  <- seq(w.hat - dw, w.hat + dw, length=31)

## Sanity and time check
schedule1 <- schedule.from.times(res$times)
system.time(delta <- run(w.hat, p, schedule1) - w.hat)
delta # should be very close to zero

seed_rain.in  <- seq(w.hat - dw, w.hat + dw, length=31)
seed_rain.out <- unlist(mclapply(seed_rain.in, run, p, schedule1))

fit <- lm(seed_rain.out ~ seed_rain.in)

## Here is input seeds vs output seeds:
plot(seed_rain.in, seed_rain.out, xlab="Incoming seed rain",
     ylab="Outgoing seed rain", las=1)
abline(0, 1, lty=2, col="grey")
abline(fit, lty=2)
cobweb(approach)

## # 2: Compare with the reference model
get.seed_rain.out.reference <- function(path)
  unname(load.reference.output(path)$seed_rain_out)

equilibrium.seed.rain.reference <- function(p, nsteps, path,
                                            progress=FALSE,
                                            verbose=TRUE) {
  if (!evolve.is.installed())
    install.evolve()

  history <- list()
  for (i in seq_len(nsteps)) {
    run.reference(path, p, verbose=FALSE)
    res <- list(seed_rain=cbind("in"=p$seed_rain,
                  "out"=get.seed_rain.out.reference(path)))
    history <- c(history, list(res))
    seed_rain <- res[["seed_rain"]]
    change <- seed_rain[,"out"] - seed_rain[,"in"]

    p$seed_rain <- seed_rain[,"out"]
    if (verbose)
      message(sprintf("*** %d: %2.5f -> %2.5f (delta = %2.5f)", i,
                      seed_rain[,"in"], seed_rain[,"out"], change))
  }

  if (progress)
    attr(res, "progress") <- history
  res
}

run.reference <- function(w, p, path) {
  p$seed_rain <- w
  tree::run.reference(path, p, verbose=FALSE)
  get.seed_rain.out.reference(path)
}

p$seed_rain <- 1.1
path.r <- "ref-equilibrium"
res.r <- equilibrium.seed.rain.reference(p, 10, path.r, progress=TRUE)

w.hat.r <- unname(res.r[["seed_rain"]][,"out"])
approach.r <- t(sapply(attr(res.r, "progress"), "[[", "seed_rain"))

cols <- c(t="black", r="red")

r <- range(approach, approach.r)
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach,   pch=19, cex=.5, type="o", col=cols[["t"]])
cobweb(approach.r, pch=19, cex=.5, type="o", col=cols[["r"]])

m <- (w.hat + w.hat.r)/2
d <- abs(w.hat - w.hat.r)
r <- m + c(-1, 1) * d * 1.1
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach,   pch=19, cex=.5, type="o", col=cols[["t"]])
cobweb(approach.r, pch=19, cex=.5, type="o", col=cols[["r"]])


## Note that this *cannot* be run with mclapply because it's all done
## through filesystem access and system calls.
seed_rain.in.r <- seq(w.hat.r - dw, w.hat.r + dw,
                      length=length(seed_rain.in))
seed_rain.out.r <- sapply(seed_rain.in.r, run.reference, p, path.r)

fit.r <- lm(seed_rain.out.r ~ seed_rain.in.r)

plot(seed_rain.in.r, seed_rain.out.r)
abline(0, 1)
abline(fit.r, lty=2)

## Compared against tree:
matplot(cbind(seed_rain.in,  seed_rain.in.r),
        cbind(seed_rain.out, seed_rain.out.r),
        xlab="Incoming seed rain", ylab="Outgoing seed rain",
        las=1, pch=1, col=cols)
abline(0, 1)
abline(fit,   lty=2, col=cols[["t"]])
abline(fit.r, lty=2, col=cols[["r"]])

matplot(cbind(seed_rain.in, seed_rain.in.r),
        cbind(resid(fit),   resid(fit.r)),
        xlab="Incoming seed rain", ylab="Residual seed rain",
        las=1, pch=1, col=cols)
abline(h=0, v=c(w.hat, w.hat.r))

## TODO: This is a concern -- and different to what I was seeing
## before, I think.  Cohort refinement is off, but something is up.
matplot(cbind(seed_rain.in, seed_rain.in.r),
        cbind(resid(fit),   resid(fit.r)),
        xlab="Incoming seed rain", ylab="Residual seed rain",
        las=1, pch=1, col=cols, ylim=range(resid(fit)))
abline(h=0, v=c(w.hat, w.hat.r))

## # Global function shape
seed_rain.in.global <- seq(1, w.hat + 10, length=51)

seed_rain.out.global <-
  unlist(mclapply(seed_rain.in.global, run.new.schedule, p, schedule0))
seed_rain.out.global.a <-
  unlist(lapply(seed_rain.in.global, run.reference, p, path.r))

matplot(seed_rain.in.global,
        cbind(seed_rain.out.global, seed_rain.out.global.r),
        las=1, type="l", col=cols, lty=1,
        xlab="Incoming seed rain", ylab="Outgoing seed rain")
abline(0, 1, lty=2, col="grey")
abline(fit,   lty=2, col=cols[["t"]])
abline(fit.a, lty=2, col=cols[["a"]])
cobweb(approach,   col=cols[["t"]])
cobweb(approach.a, col=cols[["a"]])

## TODO: Review the list of control parameters and see which are
## causing the problems.

## TODO: There is still some major-ish instability in the lhs of the
## global plot.  The might give some clues about where the remaining
## problems are.

## TODO: See how cohort refinement interacts with the code above
## (i.e. run an additional series where cohorts are refined at the
## same time, perhaps from the previous starting point).

## TODO: Now that output is more stable look at cohort refinement
## again

## TODO: Cohort merging, still.
