## TODO: approximate assimilation is not working due to spline
## refining failure.

## TODO: Think about if we converge to the same place with both
## approaches!

## TODO: Cohort merging needs implementing (?)

## TODO: See how this now behaves without adaptive integration.  It's
## possible that we can force this down below 1e-3 fairly easily.

library(tree)
library(parallel)

## Same set up as `splitting_cohorts.R`:
p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- 1.1                       # Whatever
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

t.linear <- seq(0, 104, length=31)
t.default <- cohort.introduction.times(104)

times.linear <- build.schedule(p, t.linear, 20, 1e-3,
                               progress=TRUE, verbose=TRUE)
times.default <- build.schedule(p, t.default, 20, 1e-3,
                                progress=TRUE, verbose=TRUE)

## Next, look at what an "optimal" schedule would look like from the
## point of only fitness looks like.  We'll aim for a schedule of 250
## total points (this is about 32 hours work, embarassingly enough).
n.total <- 250
n.linear <- n.total - length(t.linear)
n.default <- n.total - length(t.default)

n.linear <- n.default <- 10
times.w.linear <-
  run.cached(build.schedule.fitness(p, t.linear, n.linear,
                                    progress=TRUE, verbose=TRUE),
             "times.w.linear.rds")
times.w.default <-
  run.cached(build.schedule.fitness(p, t.default, n.default,
                                    progress=TRUE, verbose=TRUE),
             "times.w.default.rds")

## These look like they're grinding away and constantly increasing
## fitness.  I wonder if this is because of errors propagating up from
## lower level calculations?
w.linear  <- sapply(attr(times.w.linear,  "progress")[-1], "[[", "w")
w.default <- sapply(attr(times.w.default, "progress")[-1], "[[", "w")

w.linear.cmp <- w.linear[seq(to=length(w.linear), length=length(w.default))]

matplot(cbind(w.default, w.linear.cmp),
        type="o", col=c("red", "blue"), pch=19, cex=.3, lty=1)

matplot(cbind(w.default, w.linear.cmp),
        ylim=c(max(w.default)-1, max(w.default)),
        type="o", col=c("red", "blue"), pch=19, cex=.3, lty=1)

rm(p)

## For two species:
p2 <- new(Parameters)
p2$add_strategy(new(Strategy, list(lma=0.0648406, hmat=26.3098)))
p2$add_strategy(new(Strategy, list(lma=0.1977910, hmat=27.8790)))
p2$seed_rain <- c(1.1, 2.1)               # Whatever
p2$set_parameters(list(patch_area=1.0))   # See issue #13
p2$set_control_parameters(fast.control()) # A bit faster

t.linear <- seq(0, 104, length=31)
t.default <- cohort.introduction.times(104)

times.linear <- build.schedule(p2, t.linear, 20, 1e-3,
                               progress=TRUE, verbose=TRUE)
times.default <- build.schedule(p2, t.default, 20, 1e-3,
                                progress=TRUE, verbose=TRUE)
