library(tree)

## Some functions for exploring introduction times.
##
## This is similar to tree::cohort_introduction_times(), but does not
## do the `2^floor(log2())` transformation.
cohort_introduction_times <- function(max.time, multiplier=0.2,
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

insert_time <- function(i, x) {
  j <- seq_len(i)
  c(x[j], (x[i] + x[i+1])/2, x[-j])
}

run_with_times <- function(times, ebt) {
  ebt$reset()
  ebt$set_times(times, 1L)
  ebt$run()
  ebt$fitness(1L)
}

## This should really move into tree
cohort_fitness <- function(ebt) {
  cbind(t=ebt$cohort_schedule$times(1),
        seeds=ebt$fitness_cohort(1))
}

p <- new(Parameters)
p$add_strategy(new(Strategy))

p$seed_rain <- 1.1
p$set_parameters(list(patch_area=1.0)) # See issue #13

p$set_control_parameters(fast_control())

## # 1. The fitness calculation depends on the cohort spacing

## Three progressively more closely spaced times:
eps <- p$control$parameters$schedule_default_patch_survival
t_max <- p$disturbance$cdf(eps)
tt_1 <- cohort_introduction_times(t_max)
tt_2 <- interleave(tt_1)
tt_3 <- interleave(tt_2)

ebt <- new(EBT, p)
ebt$cohort_schedule$max_time <- t_max

w_1 <- run_with_times(tt_1, ebt)
w_2 <- run_with_times(tt_2, ebt)
w_3 <- run_with_times(tt_3, ebt)

## Fitness increases as time is cohorts are introduced more finely,
## though at a potentially saturating rate.  We're doing lots more
## work at the more refined end though!
c(w_1, w_2, w_3)

## # 2: Where is this difference coming from?

## Are we just geting more accurate resolution of some of the time
## courses?  Is it mostly coming from the beginning individuals or the
## end ones?

## By brute force, split cohorts along the time course and see how we
## do; that might give some clue:
run_with_insert <- function(i, t, ebt) {
  run_with_times(insert_time(i, t), ebt)
}

## This takes a little while because we're re-running the entire
## simulation (about 3 minutes).  This would be a bit easier to do if
## we could easily restart the entire simulation because we could just
## to that for the point just before the new introduction.  See issue
## #40.
##+ cache=TRUE
i <- seq_len(length(tt_1) - 1)
res <- sapply(i, run_with_insert, tt_1, ebt)

## This is really interesting; refining most cohorts does not change
## anything, but there is a little period around the 47th cohort where
## it makes a massive difference, increasing fitness by about 30,
## which is most of the difference that we saw between the refined and
## non-refined sets (total increase from doubling resolution is 52)
dw_1 <- res - w_1
tm_1 <- (tt_1[-1] + tt_1[-length(tt_1)])/2

##+ fitness_difference_by_index
plot(dw_1, xlab="Cohort insertion point", ylab="Fitness difference")
##+ fitness_difference_by_time
plot(tm_1, dw_1, xlab="Cohort insertion time",
     ylab="Fitness difference")
##+ fitness_difference_by_time_log
plot(tm_1, dw_1, xlab="Cohort insertion time",
     ylab="Fitness difference", log="x")

## The reason for the change in total fitness is how we compute the
## contribution to fitness.  Pulling apart the fitness calculations
## above to get the per-cohort contributions:
tmp <- cohort_fitness(ebt)

## And plotting this with the position (vertical line) of the most
## influential cohort introduction time.  Dashed vertical lines are
## the 5 next most important splits, from darkest (2nd most important)
## to lightest (6th most important).
##+ fitness_difference_most_influential
plot(tmp[-1,], log="x", las=1,
     xlab="Introduction time", ylab="Fitness contribution")
abline(v=tm_1[order(abs(dw_1), decreasing=TRUE)[2:6]],
       col=grey.colors(5), lty=2)
abline(v=tm_1[which.max(abs(dw_1))], col="red")

## Characterising the decrease in fitness seems to be the important
## step, especially as the derivative of fitness with time increases
## and as

## Here is the same figure on a non-log basis, covering the important
## range of times.
##+ fitness_difference_most_influential_non_log
r <- range(tm_1[order(abs(dw_1), decreasing=TRUE)[1:6]])
plot(tmp[-1,], las=1, xlim=r,
     xlab="Introduction time", ylab="Fitness contribution")
abline(v=tm_1[order(abs(dw_1), decreasing=TRUE)[2:6]],
       col=grey.colors(5), lty=2)
abline(v=tm_1[which.max(abs(dw_1))], col="red")

## 3: What is it about this point that is important?
idx <- which.max(abs(dw_1))

sched <- new(CohortSchedule, 1)
sched$max_time <- max(tt_1)
sched$set_times(tt_1, 1)

res_0 <- run_ebt_collect(ebt$parameters, sched)

sched$set_times(insert_time(idx, tt_1), 1)
res_cmp <- run_ebt_collect(ebt$parameters, sched)

## Confirm where the new cohort is (at idx)
identical(res_0$times, res_cmp$times[-idx])

## Look at the light environments over time
env_0   <- res_0$light_env
env_cmp <- res_cmp$light_env[-idx]

show_env <- function(i) {
  plot(env_0[[i]], type="o")
  lines(env_cmp[[i]], col="red", type="o", cex=.5, lty=2)
}

## These give basically the same light environment until 74 -- at this
## point they even have the same dimensions:
n <- cbind(sapply(env_0, nrow),
           sapply(env_cmp, nrow))
idx2 <- which(n[,1] != n[,2])[1]

## The light environments are basically the same looking the whole way
## through here, though refined quite differently:

##+ env_5
show_env(idx2 + 5)
##+ env_10
show_env(idx2 + 10)
##+ env_15
show_env(idx2 + 15)
## that difference persists:
##+ env_20
show_env(idx2 + 20)
##+ env_30
show_env(idx2 + 30)
##+ env_40
show_env(idx2 + 40)
##+ env_final
show_env(length(env_0))
