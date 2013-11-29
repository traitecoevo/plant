## TODO: approximate assimilation is not working due to spline
## refining failure.

## TODO: Think about what function of dt/t looks like (it is
## approximately log-linear, excluding the first cohorts) and work out
## how to quickly jump approximately to the right shape.  Dig out the
## heuristic from Daniel's code a second time.

## TODO: Think about if we converge to the same place with both
## approaches!

## TODO: Compare with Daniel's approach?

## TODO: Cohort merging needs implementing (?)
library(tree)
library(parallel)

source("build_schedule-fun.R")

## Same set up as `splitting_cohorts.R`:
p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- 1.1                       # Whatever
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

times <- build.schedule(p, 20, 30, 104, 1e-3,
                        progress=TRUE, verbose=TRUE)

## Look at the error at the last step:
res <- attr(times, "progress")
suppressWarnings(plot(last(res)$err$total, type="o", log="y"))
abline(h=1e-3)

## Next, look at what an "optimal" schedule would look like from the
## point of only fitness looks like.
times.w <- run.cached(build.schedule.fitness(p, 200, 30, 104, 4,
                                             progress=TRUE,
                                             verbose=TRUE),
                      "schedule.fitness.rds")
