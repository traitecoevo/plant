library(tree)
library(parallel)

## Try to establish what the equilubrium seed rain is.  We won't
## actually do this very often, because we'll rely on stochastic
## assembly to do this for us.
p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- 1.1                       # Starting rain.
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

## TODO: Still have the issue that the times created here don't appear
##   to reach 104 (that might actually account for some of the
##   difference with the reference model?)
times <- build.schedule(p, 20, cohort.introduction.times(104), 1e-3,
                        progress=TRUE, verbose=TRUE)
res <- list(seeds.in  = p$seed_rain,
            seeds.out = last(attr(times, "progress"))$w,
            times     = as.numeric(times))

history <- list(res)

approach <- t(sapply(history, function(x)
                     x[c("seeds.in", "seeds.out")]))
colnames(approach) <- c("in", "out")

## Once we have the cohort merging working, we could use the previous
## iteration's schedule here, I think.  Plus do a savage cull before
## starting?

## Iterate the seed rain in/out until we somewhat stabilise.  Oddly
## with the switch to a
for (i in seq_len(5)) {
  p$seed_rain <- last(history)$seeds.out
  times <- build.schedule(p, 20, cohort.introduction.times(104), 1e-3,
                          progress=TRUE, verbose=TRUE)
  res <- list(seeds.in  = p$seed_rain,
              seeds.out = last(attr(times, "progress"))$w,
              times     = as.numeric(times))
  history <- c(history, list(res))
  message(sprintf("*** %d: %2.5f (%2.5f)", i,
                  res$seeds.out, res$seeds.out - res$seeds.in))
}

## Then, in the vinicity of the root we should look at what the curve
## actually looks like, without adaptive refinement.
f <- function(w, p, times) {
  p$seed_rain <- w
  ebt <- new(EBT, p)
  ebt$cohort_schedule <- schedule.from.times(times)
  ebt$run()
  ebt$fitness(1)
}

w.hat <- last(history)[["seeds.out"]]
times <- last(history)[["times"]]
dw <- 2 # range of input to vary (plus and minus this many seeds)

seeds.in <- seq(w.hat - dw, w.hat + dw, length=31)
seeds.out <- unlist(mclapply(seeds.in, f, p, times))

fit <- lm(seeds.out ~ seeds.in)

## Here is input seeds vs output seeds; this function *should* be
## smooth.
plot(seeds.in, seeds.out, xlab="Incoming seed rain",
     ylab="Outgoing seed rain", las=1)
abline(0, 1)
abline(fit, lty=2)

## TODO: This shows that we're off by a bit because of not including the
## last time.  Gah.  Need to go through and sort times out soon.
##
##   points(approach)
##
## TODO: Once that's worked out, add the approach points on here.

########################################

## Rerun, but allowing adaptive integration during the calculation:
p.adaptive <- p$clone()
p.adaptive$set_control_parameters(list(plant_assimilation_adaptive=TRUE))

## This runs quite a bit slower:
system.time(f(w.hat, p,          times)) #  9.1s
system.time(f(w.hat, p.adaptive, times)) # 18.5s

seeds.out.adaptive <- unlist(mclapply(seeds.in, f, p.adaptive, times))
fit.adaptive <- lm(seeds.out.adaptive ~ seeds.in)

matplot(cbind(seeds.in,  seeds.in),
        cbind(seeds.out, seeds.out.adaptive),
        xlab="Incoming seed rain", ylab="Outgoing seed rain",
        las=1, pch=1, col=c("black", "red"))
abline(0, 1)
abline(fit,          lty=2)
abline(fit.adaptive, lty=2, col="red")

## This might be easier to see in terms of residuals:

matplot(seeds.in,
        cbind(resid(fit), resid(fit.adaptive)),
        xlab="Incoming seed rain", ylab="Residual seed rain",
        las=1, pch=1, col=c("black", "red"))
abline(h=0, v=w.hat)

########################################

## Rerun, using a more fine-scale error control, especially around the
## assimilation calculation, but keeping the non-adaptive
## integration.
p.fine <- p$clone()
p.fine$set_control_parameters(list(plant_assimilation_tol=1e-6,
                                   ode_tol_rel=1e-6,
                                   ode_tol_abs=1e-6))

## This runs quite a bit slower:
system.time(f(w.hat, p,      times)) #  9.1s
system.time(f(w.hat, p.fine, times)) # 46.0s

seeds.out.fine <- unlist(mclapply(seeds.in, f, p.fine, times))

matplot(seeds.in,
        cbind(seeds.out, seeds.out.adaptive, seeds.out.fine),
        xlab="Incoming seed rain", ylab="Outgoing seed rain",
        las=1, pch=1, col=c("black", "red", "blue"))
abline(0, 1)
abline(fit,          lty=2)
abline(fit.adaptive, lty=2, col="red")
abline(fit.fine,     lty=2, col="blue")

fit.fine     <- lm(seeds.out.fine     ~ seeds.in)
matplot(seeds.in,
        cbind(resid(fit), resid(fit.adaptive), resid(fit.fine)),
        xlab="Incoming seed rain", ylab="Residual seed rain",
        las=1, pch=1, col=c("black", "red", "blue"))
abline(h=0, v=w.hat)

########################################

## TODO: Rerun, but allowing the cohort schedule to be built during
## the approach.

########################################

## Compare with falster-traitdiversity.  This is a bit more
## complicated because the cohort refinement is going on as well
## (whereas above we are using the same spacings).  So we don't really
## have comparable levels of "noise".
if (!evolve.is.installed())
  install.evolve()

get.fitness.reference <- function(path)
  unname(load.reference.output(path)$seed_rain_out)

## Because of slight differences, we don't equilibrate in the same
## place, so need to home in on the right spot again:
path <- "ref-equilibrium"
p$seed_rain <- 1.1
run.reference(path, p, verbose=FALSE)
res <- list(seeds.in  = p$seed_rain,
            seeds.out = get.fitness.reference(path))

history.ref <- list(res)
for (i in seq_len(5)) {
  p$seed_rain <- last(history.ref)$seeds.out
  run.reference(path, p, verbose=FALSE)
  res <- list(seeds.in  = p$seed_rain,
              seeds.out = get.fitness.reference(path))
  history.ref <- c(history.ref, list(res))
  message(sprintf("*** %d: %2.5f (%2.5f)", i,
                  res$seeds.out, res$seeds.out - res$seeds.in))
}

w.hat.ref <- last(history.ref)[["seeds.out"]]
f.ref <- function(w, p) {
  p$seed_rain <- w
  run.reference(path, p, verbose=FALSE)
  get.fitness.reference(path)
}

## Note that this *cannot* be run with mclapply because it's all done
## through filesystem access and system calls.
seeds.in.ref <- seq(w.hat.ref - dw, w.hat.ref + dw, length=31)
seeds.out.ref <- sapply(seeds.in.ref, f.ref, p)

fit.ref <- lm(seeds.out.ref ~ seeds.in.ref)

plot(seeds.in.ref, seeds.out.ref)
abline(0, 1)
abline(fit.ref, lty=2)

## Compared against tree:
matplot(cbind(seeds.in,  seeds.in,           seeds.in.ref),
        cbind(seeds.out, seeds.out.adaptive, seeds.out.ref),
        xlab="Incoming seed rain", ylab="Outgoing seed rain",
        las=1, pch=1, col=c("black", "red", "green4"))
abline(0, 1)
abline(fit,          lty=2)
abline(fit.adaptive, lty=2, col="red")
abline(fit.ref,      lty=2, col="green4")

matplot(cbind(seeds.in,   seeds.in,            seeds.in.ref),
        cbind(resid(fit), resid(fit.adaptive), resid(fit.ref)),
        xlab="Incoming seed rain", ylab="Residual seed rain",
        las=1, pch=1, col=c("black", "red", "green4"))
abline(h=0, v=w.hat)

########################################

## Next, run a far more extensive range of incoming seed rain:
seeds.in.wide <- seq(1, max(w.hat, w.hat.ref) + 10, length=51)
## Without changing the cohort spacing (might not be very clever)
## (also note that these are going to be run with the "fine" times,
## not the base times).
seeds.out.wide <- unlist(mclapply(seeds.in.wide, f, p, times))
seeds.out.wide.ref <- sapply(seeds.in.wide, f.ref, p)

matplot(seeds.in.wide, cbind(seeds.out.wide, seeds.out.wide.ref),
        las=1, type="l", col=c("black", "blue"), lty=1,
        xlab="Incoming seed rain", ylab="Outgoing seed rain")
abline(0, 1, col="darkgrey")
abline(fit,          lty=2)
abline(fit.adaptive, lty=2, col="red")
abline(fit.ref,      lty=2, col="green4")

matplot(seeds.in.wide, cbind(seeds.out.wide, seeds.out.wide.ref),
        las=1, type="l", col=c("black", "blue"), lty=1, asp=1,
        xlab="Incoming seed rain", ylab="Outgoing seed rain")
abline(0, 1, col="darkgrey")
abline(fit,          lty=2)
abline(fit.adaptive, lty=2, col="red")
abline(fit.ref,      lty=2, col="green4")

## TODO: Review the list of control parameters and see which are
## causing the problems.

## TODO: There is still some major-ish instability in the lhs of the
## graph.  The might give some clues about where the remaining
## problems are.

## TODO: See how cohort refinement interacts with the code above
## (i.e. run an additional series where cohorts are refined at the
## same time, perhaps from the previous starting point).

## TODO: Now that output is more stable look at cohort refinement
## again

## TODO: Cohort merging, still.
