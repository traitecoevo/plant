source("helper-tree.R")

context("EBTMutantRunner")

p <- new(Parameters)
p$add_strategy(new(Strategy))
p$set_control_parameters(fast.control())
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$seed_rain <- 1.1
t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)

## This takes a few seconds to run
sched <- schedule.from.times(cohort.introduction.times(t.max))
res <- run.ebt.collect(p, sched)

ebt <- new(EBT, p)
ebt$cohort_schedule <- sched
system.time(ebt$run()) # 0.716

f <- function(xy) {
  s <- new(Interpolator)
  s$init(xy[,1], xy[,2])
  s
}

## Almost the entire time here is spent building the set of splines;
## they were already in one piece before ebt.run.collect pulled them
## apart (perhaps it should leave them together?).  This is not used
## so often.
obj <- new(FakeLightEnvironment, res$t, lapply(res$light.env, f))

p.m <- p$copy()
expect_that(new(tree:::EBTMutantRunner, p.m, obj),
            throws_error())
p.m$set_control_parameters(list(environment_light_skip=TRUE))

ebt.m <- new(tree:::EBTMutantRunner, p.m, obj)

g <- function(ebt)
  ebt$patch$environment$light_environment

## First, try with the same schedule.
ebt.m$cohort_schedule <- sched

system.time(ebt.m$run()) # 0.415

## Surprisingly different:
ebt.m$fitnesses # 636.4227
ebt$fitnesses   # 563.8381

quad <- new(QK, 51)
t.int <- quad$integrate_vector_x(0, t.max)

tt <- c(0, sort(t.int), t.max)
sched.m <- sched$copy()
sched.m$set_times(tt, 1)

ebt.m <- new(tree:::EBTMutantRunner, p.m, obj)
ebt.m$cohort_schedule <- sched.m
ebt.m$run() # 0.124 s

y <- ebt.m$fitness_cohort(1)

## Now, just get y back into the right order and we should be able to
## feed these back into the integrator.

yy <- y[-c(1, length(y))]
i <- seq_along(t.int)[order(t.int)]
w <- quad$integrate_vector(yy[i], 0, t.max) # 728.6976
ebt.m$fitnesses # 968.7117
