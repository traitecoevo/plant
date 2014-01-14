library(tree)

## Some functions for exploring introduction times.
##
## This is similar to tree::cohort.introduction.times(), but does not
## do the `2^floor(log2())` transformation.
cohort.introduction.times <- function(max.time, multiplier=0.2,
                                      min.step.size=1e-05,
                                      max.step.size=2) {
  if (min.step.size <= 0)
    stop("The minimum step size must be greater than zero")
  trim <- function(x)
    max(min(x, max.step.size), min.step.size)
  times <- numeric(0)
  dt <- time <- 0
  while (time <= max.time) {
    times <- c(times, time)
    time <- time + trim(time * multiplier)
  }
  times
}

interleave <- function(x) {
  n <- length(x)
  xp <- c(x[-n] + diff(x)/2, NA)
  c(rbind(x, xp))[-2*n]
}

insert.time <- function(i, x) {
  j <- seq_len(i)
  c(x[j], (x[i] + x[i+1])/2, x[-j])
}

run.with.times <- function(times, ebt) {
  ebt$reset()
  ebt$set_times(times, 1L)
  ebt$run()
  ebt$fitness(1L)
}

## This should really move into tree
cohort.fitness <- function(ebt) {
  cbind(t=ebt$cohort_schedule$times(1),
        seeds=ebt$fitness_cohort(1))
}

p <- new(Parameters)
p$add_strategy(new(Strategy))

p$seed_rain <- 1.1
p$set_parameters(list(patch_area=1.0)) # See issue #13

## Relatively quick control settings:
ctrl.new <- list()
ctrl.new$environment_light_rescale_usually <- TRUE
ctrl.new$environment_light_tol <- 1e-4
ctrl.new$plant_assimilation_rule <- 21
ctrl.new$plant_assimilation_over_distribution <- FALSE
ctrl.new$plant_assimilation_tol <- 1e-4
ctrl.new$ode_tol_rel <- 1e-4
ctrl.new$ode_tol_abs <- 1e-4
ctrl.new$ode_step_size_max <- 5
ctrl.new$cohort_gradient_direction <- -1
ctrl.new$cohort_gradient_richardson <- FALSE
p$set_control_parameters(ctrl.new)

ebt <- new(EBT, p)

## # 1. The fitness calculation depends on the cohort spacing

## Three progressively more closely spaced times:
t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
tt.1 <- cohort.introduction.times(t.max)
tt.2 <- interleave(tt.1)
tt.3 <- interleave(tt.2)

w.1 <- run.with.times(tt.1, ebt)
w.2 <- run.with.times(tt.2, ebt)
w.3 <- run.with.times(tt.3, ebt)

## Fitness increases as time is cohorts are introduced more finely,
## though at a potentially saturating rate.  We're doing lots more
## work at the more refined end though!
c(w.1, w.2, w.3)

## # 2: Where is this difference coming from?

## Are we just geting more accurate resolution of some of the time
## courses?  Is it mostly coming from the beginning individuals or the
## end ones?

## By brute force, split cohorts along the time course and see how we
## do; that might give some clue:
run.with.insert <- function(i, t, ebt)
  run.with.times(insert.time(i, t), ebt)

## This takes a little while because we're re-running the entire
## simulation (about 3 minutes).  This would be a bit easier to do if
## we could easily restart the entire simulation because we could just
## to that for the point just before the new introduction.  See issue
## #40.
i <- seq_len(length(tt.1) - 1)
res <- sapply(i, run.with.insert, tt.1, ebt)

## This is really interesting; refining most cohorts does not change
## anything, but there is a little period around the 47th cohort where
## it makes a massive difference, increasing fitness by about 30,
## which is most of the difference that we saw between the refined and
## non-refined sets (total increase from doubling resolution is 52)
dw.1 <- res - w.1
tm.1 <- (tt.1[-1] + tt.1[-length(tt.1)])/2

plot(dw.1, xlab="Cohort insertion point", ylab="Fitness difference")
plot(tm.1, dw.1, xlab="Cohort insertion time",
     ylab="Fitness difference")
plot(tm.1, dw.1, xlab="Cohort insertion time",
     ylab="Fitness difference", log="x")

## The reason for the change in total fitness is how we compute the
## contribution to fitness.  Pulling apart the fitness calculations
## above to get the per-cohort contributions:
tmp <- cohort.fitness(ebt)

## And plotting this with the position (vertical line) of the most
## influential cohort introduction time.  Dashed vertical lines are
## the 5 next most important splits, from darkest (2nd most important)
## to lightest (6th most important).
plot(tmp[-1,], log="x", las=1,
     xlab="Introduction time", ylab="Fitness contribution")
abline(v=tm.1[order(abs(dw.1), decreasing=TRUE)[2:6]],
       col=grey.colors(5), lty=2)
abline(v=tm.1[which.max(abs(dw.1))], col="red")

## Characterising the decrease in fitness seems to be the important
## step, especially as the derivative of fitness with time increases
## and as

## Here is the same figure on a non-log basis, covering the important
## range of times.
r <- range(tm.1[order(abs(dw.1), decreasing=TRUE)[1:6]])
plot(tmp[-1,], las=1, xlim=r,
     xlab="Introduction time", ylab="Fitness contribution")
abline(v=tm.1[order(abs(dw.1), decreasing=TRUE)[2:6]],
       col=grey.colors(5), lty=2)
abline(v=tm.1[which.max(abs(dw.1))], col="red")

## 3: What is it about this point that is important?
idx <- which.max(abs(dw.1))

sched <- new(CohortSchedule, 1)
sched$max_time <- max(tt.1)
sched$set_times(tt.1, 1)

res.0 <- run.ebt.collect(ebt$parameters, sched)

sched$set_times(insert.time(idx, tt.1), 1)
res.cmp <- run.ebt.collect(ebt$parameters, sched)

## Confirm where the new cohort is (at idx)
identical(res.0$times, res.cmp$times[-idx])

## Look at the light environments over time
env.0   <- res.0$light.env
env.cmp <- res.cmp$light.env[-idx]

show.env <- function(i) {
  plot(env.0[[i]], type="o")
  lines(env.cmp[[i]], col="red", type="o", cex=.5, lty=2)
}

## These give basically the same light environment until 74 -- at this
## point they even have the same dimensions:
n <- cbind(sapply(env.0, nrow),
           sapply(env.cmp, nrow))
idx2 <- which(n[,1] != n[,2])[1]

## Differences by this point are fairly small
all.equal(env.0[idx:(idx2 - 1)],
          env.cmp[idx:(idx2 - 1)], tolerance=1e-4)

## After this point, the version with the additional cohort has
## *greater* canopy openness at the bottom of the light environment
## (so below 10m tall, there is more light getting through).
show.env(idx2 + 5)

## This becomes a wedge, with another possible round of recruitment,
## causing a second drop in light at about 3m:
show.env(idx2 + 9)

## By the time a similar wedge starts appearing in the first version,
## the wedge in the new verion has steepend dramatically -- now the
## light environment is lower for very small heights than for medium
## heights:
show.env(idx2 + 12)

## Then, the light environment at the forest floor equalises for the
## two versions, but the position of the wedge varies between runs:
show.env(idx2 + 15)

## that difference persists:
show.env(idx2 + 20)
show.env(idx2 + 30)
show.env(idx2 + 40)

## The final light environment is fairly similar though
show.env(length(env.0))
