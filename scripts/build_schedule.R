library(tree)
library(parallel)

## Same set up as `splitting_cohorts.R`:
p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- 1.1                       # Whatever
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast_control()) # A bit faster

## Useful for now, given the idea is to see how schedules are built:
p$set_control_parameters(list(schedule_verbose=TRUE))

## Build a schedule based on the parameters:
sched_default <- build_schedule(p)

## Again, with a linear initial schedule:
sched_linear0 <- new(CohortSchedule, 1L)
sched_linear0$max_time <- sched_default$max_time
sched_linear0$set_times(seq(0, sched_linear0$max_time, length=31), 1L)
sched_linear <- build_schedule(p, sched_linear0)

## For two species:
p2 <- new(Parameters)
p2$add_strategy(new(Strategy, list(lma=0.0648406, hmat=26.3098)))
p2$add_strategy(new(Strategy, list(lma=0.1977910, hmat=27.8790)))
p2$seed_rain <- c(1.1, 2.1)               # Whatever
p2$set_parameters(list(patch_area=1.0))   # See issue #13
p2$set_control_parameters(fast_control()) # A bit faster
p2$set_control_parameters(list(schedule_verbose=TRUE))

sched2_default <- build_schedule(p2)
