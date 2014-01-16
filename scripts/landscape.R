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
schedule0 <- schedule.from.times(cohort.introduction.times(t.max))
schedule <- build.schedule(p, schedule0, 20, 1e-2,
                           progress=FALSE, verbose=TRUE)

## Once the schedule building is improved, we'll fix it so that this
## is not required.  We should return a list with all the required
## bits I think.

## Running this takes about 2.9s.
ebt <- run.ebt(p, schedule$copy())
schedule$ode_times <- ebt$ode_times

## Now, let's get a view of the fitness landscape with respect to one
## parameter.
lma <- p[[1]]$parameters[["lma"]]

n.mutants <- 11

## Add the mutants; they will all be introduced with a seed rain of 1,
## though they will not influence the environent at all.
p.with.mutants <- p$copy()
## Add a vector of new lmas:
p.with.mutants$add_strategy_mutant(new(Strategy, list(lma=lma)))
for (i in seq(lma * 0.1, lma * 1.2, length.out=n.mutants))
  p.with.mutants$add_strategy_mutant(new(Strategy, list(lma=i)))

## Probably could benefit from an expand method?
##
## TODO: simple n_mutant and n_resident method for Parameters
expand.schedule <- function(schedule, n.mutant) {
  n.resident <- schedule$n_species
  if (n.resident != 1)
    stop("Can only expand a schedule with one resident at the moment")
  ret <- new(CohortSchedule, n.resident + n.mutant)
  ret$max_time <- schedule$max_time

  ## Copy residents over:
  for (i in seq_len(n.resident))
    ret$set_times(schedule$times(i), i)

  ## TODO: work out what to do better here:
  times.mutant <- schedule$times(1)
  for (i in seq_len(n.mutant))
    ret$set_times(times.mutant, n.resident + i)

  ret$ode_times <- schedule$ode_times

  ret
}

## Next, update the schedule!
## TODO: update to do p$n_mutants
schedule2 <- expand.schedule(schedule, n.mutants+1)

## And run, producing only fitnesses.  The resident fitness should be
## unchanged here.
##
## Unfortunately this is very slow!  We probably need to think
## carefully about how to speed this up.
ebt.with.mutants <- run.ebt(p.with.mutants, schedule2)
w.with.mutants <- ebt.with.mutants$fitnesses

lma.v <- sapply(seq_len(p.with.mutants$size),
                function(i) p.with.mutants[[i]]$parameters[["lma"]])

plot(lma.v[-(1:2)], w.with.mutants[-(1:2)], log="y")
points(lma.v[[1]], w.with.mutants[[1]] / p$seed_rain[[1]],
       col="red", pch=19)
points(lma.v[[2]], w.with.mutants[[2]],
       col="blue", pch=19, cex=.5)
