## Towards evolutionary assembly.
library(tree)
library(assertthat)
library(lhs)
library(mvtnorm)

source("assembly-fun.R")

bounds <- rbind(lma=10^c(-2.5, 1.5),
                hmat=c(0.1,100))
colnames(bounds) <- c("lower", "upper")

## VCV that explores fraction p of the bounds.
p <- 0.001
vcv <- diag(2) * p * as.numeric(diff(t(log(bounds))))

new.phenotypes <- make.new.phenotypes(2, vcv, 1, bounds)

## set.seed(1)
## traits  <- rbind(rowMeans(bounds))
## repeat {
##   plot(traits, xlim=bounds[1,], ylim=bounds[2,], log="xy")
##   traits.new <- f(traits, rep(1, nrow(traits)))
##   points(traits.new, col="blue", cex=.5, pch=19)
##   traits <- rbind(traits, traits.new)
## }

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast.control())
t.max <- p0$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)

rain.intro <- 1e-3 # rain for new arrivals
rain.eps   <- 1e-4 # below which, extinct

traits <- cbind(lma=0.1978791, hmat=16.5958691)
rain.in <- 505.55407
p <- setup.parameters(traits, rain.in, p0)

## Build a schedule for the first trait:
schedule0 <- schedule.from.times(times0, 1)
schedule <- build.schedule(p, schedule0, 20, 1e-3, progress=FALSE,
                           verbose=TRUE)

times <- schedule$all_times
rain.out <- unname(attr(schedule, "seed_rain", exact=TRUE)[,"out"])

## Still have to drop the things that don't make it at this step.
set.seed(1)
traits.new  <- new.phenotypes(traits, rain.out)
n.new       <- nrow(traits.new)
rain.in.new <- rep(rain.intro, n.new)
times.new   <- rep(list(times0), n.new)

p2 <- setup.parameters(rbind(traits, traits.new), c(rain.in, rain.in.new),
                       p0)
schedule2 <- setup.schedule(c(times, times.new), t.max)

## This is really going bonkers on cohort refinement -- way more than
## I'd have thought given that the individuals are coming in at such a
## low frequency.  The two traits that we did put in are basically the
## same as the focal individual but they get refined down to
## incredible levels (almost 600 individuals).
##
## The extra times are turning up in a weird pulse around 20.
schedule2.2 <- build.schedule(p2, schedule2, 20, 1e-3, progress=FALSE,
                              verbose=TRUE)

times2 <- schedule2.2$all_times
## The fact that seed rain has been driven down to 1e-7 really
## suggests that the cohort refinement code just needs adjusting.  We
## simply do not need the precision when things are that bad.  At the
## same time it suggests why the refinement code has bailed.
rain.out2 <- unname(attr(schedule2.2, "seed_rain", exact=TRUE)[,"out"])
