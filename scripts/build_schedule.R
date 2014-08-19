library(tree)
library(parallel)

p <- ebt_base_parameters()
p$add_strategy(new(Strategy))
p$seed_rain <- 1.1

## Build a schedule based on the parameters:
sched_default <- build_schedule(p)

## Again, with a linear initial schedule:
sched_linear0 <- new(CohortSchedule, 1L)
sched_linear0$max_time <- sched_default$max_time
sched_linear0$set_times(seq(0, sched_linear0$max_time, length=31), 1L)
sched_linear <- build_schedule(p, sched_linear0)

p2 <- ebt_base_parameters()
p2$add_strategy(new(Strategy, list(lma=0.0648406, hmat=26.3098)))
p2$add_strategy(new(Strategy, list(lma=0.1977910, hmat=27.8790)))
p2$seed_rain <- c(1.1, 2.1)               # Whatever

sched2_default <- build_schedule(p2)
