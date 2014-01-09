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

times <- build.schedule(p, 20, cohort.introduction.times(104), 1e-3,
                        progress=TRUE, verbose=TRUE)
res <- list(seeds.in  = p$seed_rain,
            seeds.out = last(attr(times, "progress"))$w,
            times     = as.numeric(times))

history <- list(res)

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
diff <- sapply(history, function(x) x$seeds.out - x$seeds.in)
i <- which.min(abs(diff))
w.hat <- history[[i]]$seeds.in
times <- history[[i]]$times

f <- function(w, p, times) {
  p$seed_rain <- w
  ebt <- new(EBT, p)
  ebt$cohort_schedule <- schedule.from.times(times)
  ebt$run()
  ebt$fitness(1)
}

dw <- 2
w.hat <- last(history)[["seeds.out"]]
seeds.in <- seq(w.hat - dw, w.hat + dw, length=31)
seeds.out <- unlist(mclapply(seeds.in, f, p, times))

## This looks way better than before; that's really great.
plot(seeds.in, seeds.out, xlab="Incoming seed rain",
     ylab="Outgoing seed rain", las=1)
abline(0, 1)
abline(lm(seeds.out ~ seeds.in), lty=2)

########################################

## Rerun, using a more fine-scale error control, especially around the
## assimilation calculation.
ctrl.fine <- fast.control()
ctrl.fine$plant_assimilation_tol <- 1e-6
ctrl.fine$ode_tol_rel <- 1e-6
ctrl.fine$ode_tol_abs <- 1e-6

p.fine <- p$clone()
p.fine$set_control_parameters(ctrl)

## This runs quite a bit slower:
system.time(f(w.hat, p,      times)) #  9.1s
system.time(f(w.hat, p.fine, times)) # 46.0s

## Run over the whole range of incoming seeds (this is very slow!)
seeds.out.fine <- unlist(mclapply(seeds.in, f, p.fine, times))

matplot(cbind(seeds.in,  seeds.in),
        cbind(seeds.out, seeds.out.fine),
        xlab="Incoming seed rain", ylab="Outgoing seed rain",
        las=1, pch=1, col=c("black", "red"))
abline(0, 1)
abline(lm(seeds.out      ~ seeds.in), lty=2)
abline(lm(seeds.out.fine ~ seeds.in), lty=2, col="red")

## This might be easier to see in terms of residuals:
fit1 <- lm(seeds.out      ~ seeds.in)
fit2 <- lm(seeds.out.fine ~ seeds.in)
matplot(seeds.in, cbind(resid(fit1), resid(fit2)),
        xlab="Incoming seed rain", ylab="Residual seed rain",
        las=1, pch=1, col=c("black", "red"))
abline(h=0, v=w.hat)

########################################

## Compare with falster-traitdiversity.  This is a bit more
## complicated because teh cohort refinement is going on as well
## (whereas above we are using the same spacings).
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

plot(seeds.in.ref, seeds.out.ref)
abline(0, 1)
abline(lm(seeds.out.ref ~ seeds.in.ref), lty=2)

## Compared against tree:
matplot(cbind(seeds.in,  seeds.in,       seeds.in.ref),
        cbind(seeds.out, seeds.out.fine, seeds.out.ref),
        xlab="Incoming seed rain", ylab="Outgoing seed rain",
        las=1, pch=1, col=c("black", "red", "blue"))
abline(0, 1)
abline(lm(seeds.out      ~ seeds.in),     lty=2)
abline(lm(seeds.out.fine ~ seeds.in),     lty=2, col="red")
abline(lm(seeds.out.ref  ~ seeds.in.ref), lty=2, col="blue")

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
abline(lm(seeds.out      ~ seeds.in),     lty=2)
abline(lm(seeds.out.fine ~ seeds.in),     lty=2, col="red")
abline(lm(seeds.out.ref  ~ seeds.in.ref), lty=2, col="blue")

## TODO: Confirm that the instability is due to the adaptive
## integration.  If so we need to be careful about how that is used in
## general.

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
