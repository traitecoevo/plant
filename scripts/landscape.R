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

## TODO: It looks to me that that max.t is ending up at 100, which is
## suboptimal, here.
max.t <- 104
times <- build.schedule(p, 20, cohort.introduction.times(max.t), 1e-2,
                        progress=FALSE, verbose=TRUE)

## Once the schedule building is improved, we'll fix it so that this
## is not required.
ebt <- run.ebt(p, schedule.from.times(times))
ode.times <- ebt$ode_times

## Now, let's get a view of the fitness landscape with respect to one
## parameter.
lma <- p[[1]]$parameters[["lma"]]

## Add the mutants; they will all be introduced with a seed rain of 1,
## though they will not influence the environent at all.
p.with.mutants <- p$clone()
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
ebt.with.mutants <- run.ebt(p.with.mutants, sched)
w.with.mutants <- ebt.with.mutants$fitnesses

lma.v <- sapply(seq_len(p.with.mutants$size),
                function(i) p.with.mutants[[i]]$parameters[["lma"]])

plot(lma.v[-(1:2)], w.with.mutants[-(1:2)], log="y")
points(lma.v[[1]], w.with.mutants[[1]] / p$seed_rain[[1]],
       col="red", pch=19)
points(lma.v[[1]], w.with.mutants[[2]],
       col="blue", pch=19, cex=.5)
