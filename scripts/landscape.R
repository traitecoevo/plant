library(tree)

## This might be useful, as we see many of the traits play out on a
## log scale.
log.seq <- function(from, to, length.out)
  exp(seq(log(from), log(to), length.out=length.out))

## Starting rain here is close to the equilibrium from
## scripts/equilibrium.R -- because the function is actually not
## terrifically smooth, it's not really the equilbrium seed rain, but
## close enough for our purposes.
p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- 505.55407
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)
times <- build.schedule(p, times0, 20, 1e-2,
                        progress=FALSE, verbose=TRUE)

## Once the schedule building is improved, we'll fix it so that this
## is not required.  We should return a list with all the required
## bits I think.

## Running this takes about 2.9s; that increases to about 6.5s if the
## adaptive integration is used (see below).
ebt <- run.ebt(p, schedule.from.times(times))
ode.times <- ebt$ode_times

p.a <- p$copy()
p.a$set_control_parameters(list(plant_assimilation_adaptive=TRUE))
system.time(ebt <- run.ebt(p.a, schedule.from.times(times)))

## Now, let's get a view of the fitness landscape with respect to one
## parameter.
lma <- p[[1]]$parameters[["lma"]]

## Add the mutants; they will all be introduced with a seed rain of 1,
## though they will not influence the environent at all.
p.with.mutants <- p$copy()
## Add the resident as a check...
p.with.mutants$add_strategy_mutant(new(Strategy, list(lma=lma)))
## ...and a vector of new lmas:
for (i in log.seq(lma * 0.1, lma * 1.2, length.out=11))
  p.with.mutants$add_strategy_mutant(new(Strategy, list(lma=i)))

## Next, update the schedule!
sched <- new(CohortSchedule, p.with.mutants$size)
for (i in seq_len(p.with.mutants$size))
  sched$set_times(times[-length(times)], i)
sched$ode_times <- ode.times

## And run, producing only fitnesses.  The resident fitness should be
## unchanged here.
##
## Unfortunately this is very slow!  We probably need to think
## carefully about how to speed this up.
ebt.with.mutants <- run.ebt(p.with.mutants, sched)
w.with.mutants <- ebt.with.mutants$fitnesses

lma.v <- sapply(seq_len(p.with.mutants$size),
                function(i) p.with.mutants[[i]]$parameters[["lma"]])

plot(lma.v[-(1:2)], w.with.mutants[-(1:2)], log="y")
points(lma.v[[1]], w.with.mutants[[1]] / p$seed_rain[[1]],
       col="red", pch=19)
points(lma.v[[1]], w.with.mutants[[2]],
       col="blue", pch=19, cex=.5)
