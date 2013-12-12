library(tree)
library(parallel)

## Try to establish what the equilubrium seed rain is.  We won't
## actually do this very often, because we'll rely on stochastic
## assembly to do this for us.
source("build_schedule-fun.R")

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

## This has the risk of the cohort refinement creating variation in
## the final fitness that will prevent this stabilising.
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

seeds.in <- seq(w.hat - 10, w.hat + 10, length=31)
seeds.out <- unlist(mclapply(seeds.in, f, p, times))

## There is quite a bit of variation here:
plot(seeds.in, seeds.out, xlab="Incoming seed rain",
     ylab="Outgoing seed rain", asp=1, las=1)
abline(0, 1)
abline(lm(seeds.out ~ seeds.in), lty=2)

## Rerun, using a more fine-scale error control, especially around the
## light environment calculation and assimilation calculation.
ctrl <- list()
ctrl$environment_light_rescale_usually <- TRUE
ctrl$environment_light_tol <- 1e-6
ctrl$plant_assimilation_over_distribution <- FALSE
ctrl$plant_assimilation_tol <- 1e-6
ctrl$ode_tol_rel <- 1e-5
ctrl$ode_tol_abs <- 1e-5
ctrl$ode_step_size_max <- 5
ctrl$cohort_gradient_direction <- -1
ctrl$cohort_gradient_richardson <- FALSE

p.fine <- new(Parameters)
p.fine$add_strategy(new(Strategy))
p.fine$set_parameters(list(patch_area=1.0)) # See issue #13
p.fine$set_control_parameters(ctrl)

## This runs quite a bit slower:
system.time(f(w.hat, p,      times)) # 16.7s
system.time(f(w.hat, p.fine, times)) # 87.5s

## Run over the whole range of incoming seeds (this is very slow!)
seeds.out.fine <- unlist(mclapply(seeds.in, f, p.fine, times))

plot(seeds.in, seeds.out, xlab="Incoming seed rain",
     ylab="Outgoing seed rain", asp=1, las=1)
points(seeds.in, seeds.out.fine, col="red", pch=19)
abline(0, 1)
abline(lm(seeds.out      ~ seeds.in), lty=2)
abline(lm(seeds.out.fine ~ seeds.in), lty=2, col="red")

## Compare this with falster-traitdiversity
path <- "ref-equilibrium"
path.evolve <- locate.evolve()
get.fitness.reference <- function(path)
  unname(load.reference.output(path)$seed_rain_out)

## Because of slight differences, we don't equilibrate in the same
## place, so need to home in on the right spot again:
p$seed_rain <- 1.1
run.reference(path, p, path.evolve=path.evolve, verbose=FALSE)
res <- list(seeds.in  = p$seed_rain,
            seeds.out = get.fitness.reference(path))
history.ref <- list(res)
for (i in seq_len(5)) {
  p$seed_rain <- last(history.ref)$seeds.out
  run.reference(path, p, path.evolve=path.evolve, verbose=FALSE)
  res <- list(seeds.in  = p$seed_rain,
              seeds.out = get.fitness.reference(path))
  history.ref <- c(history.ref, list(res))
  message(sprintf("*** %d: %2.5f (%2.5f)", i,
                  res$seeds.out, res$seeds.out - res$seeds.in))
}

diff <- sapply(history.ref, function(x) x$seeds.out - x$seeds.in)
w.hat.ref <- history.ref[[which.min(abs(diff))]]$seeds.in

f.ref <- function(w, p) {
  p$seed_rain <- w
  run.reference(path, p, path.evolve=path.evolve, verbose=FALSE)
  get.fitness.reference(path)
}

## Note that this *cannot* be run with mclapply because it's all done
## through filesystem access and system calls.
seeds.in.ref <- seq(w.hat.ref - 10, w.hat.ref + 10, length=31)
seeds.out.ref <- sapply(seeds.in.ref, f.ref, p)

plot(seeds.in.ref, seeds.out.ref, asp=1)
abline(0, 1)
abline(lm(seeds.out.ref ~ seeds.in.ref), lty=2)

## Compared against tree:
plot(seeds.in, seeds.out, xlab="Incoming seed rain",
     ylab="Outgoing seed rain", asp=1, las=1,
     xlim=range(seeds.in,  seeds.in.ref),
     ylim=range(seeds.out, seeds.out.ref))
points(seeds.in,     seeds.out.fine, col="red", pch=19)
points(seeds.in.ref, seeds.out.ref,  col="blue", pch=19)
abline(0, 1)
abline(lm(seeds.out      ~ seeds.in),     lty=2)
abline(lm(seeds.out.fine ~ seeds.in),     lty=2, col="red")
abline(lm(seeds.out.ref  ~ seeds.in.ref), lty=2, col="blue")

## Next, run a far more extensive range of incoming seed rain:
seeds.in.wide <- seq(1, max(w.hat, w.hat.ref) + 10, length=51)
## Without changing the cohort spacing (might not be very clever)
## (also note that these are going to be run with the "fine" times,
## not the base times).
seeds.out.wide <- unlist(mclapply(seeds.in.wide, f, p, times))
seeds.out.wide.ref <- sapply(seeds.in.wide, f.ref, p)

plot(seeds.in.wide, seeds.out.wide, las=1, type="l",
     xlab="Incoming seed rain", ylab="Outgoing seed rain",
     ylim=range(seeds.out.wide, seeds.out.wide.ref))
lines(seeds.in.wide, seeds.out.wide.ref,  col="blue", pch=19, cex=.5)
abline(0, 1, col="darkgrey")
abline(lm(seeds.out      ~ seeds.in),     lty=2)
abline(lm(seeds.out.fine ~ seeds.in),     lty=2, col="red")
abline(lm(seeds.out.ref  ~ seeds.in.ref), lty=2, col="blue")
