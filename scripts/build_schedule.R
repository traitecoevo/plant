## TODO: approximate assimilation is not working due to spline
## refining failure.

## TODO: Think about if we converge to the same place with both
## approaches!

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

t.linear <- seq(0, 104, length=31)
t.default <- cohort.introduction.times(104)

times.linear <- build.schedule(p, 20, t.linear, 1e-3,
                               progress=TRUE, verbose=TRUE)
times.default <- build.schedule(p, 20, t.default, 1e-3,
                                progress=TRUE, verbose=TRUE)

## Next, look at what an "optimal" schedule would look like from the
## point of only fitness looks like.  We'll aim for a schedule of 250
## total points (this is about 32 hours work, embarassingly enough).
n.total <- 250
n.linear <- n.total - length(t.linear)
n.default <- n.total - length(t.default)
times.w.linear <-
  run.cached(build.schedule.fitness(p, n.linear, t.linear, 4,
                                    progress=TRUE,
                                    verbose=TRUE),
             "times.w.linear.rds")
times.w.default <-
  run.cached(build.schedule.fitness(p, n.default, t.default, 4,
                                    progress=TRUE,
                                    verbose=TRUE),
             "times.w.default.rds")

## These look like they're grinding away and constantly increasing
## fitness.  I wonder if this is because of errors propagating up from
## lower level calculations?
w.linear  <- sapply(attr(times.w.linear,  "progress")[-1], "[[", "w")
w.default <- sapply(attr(times.w.default, "progress")[-1], "[[", "w")

plot(w.default, type="l")
lines(w.linear[seq(to=length(w.linear), length=length(w.default))])
