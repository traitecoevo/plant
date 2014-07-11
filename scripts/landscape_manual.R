## Manually construct a simple fitness landscape.  Does not actually
## use the tree::landscape function.
library(tree)

## Starting rain here is close to the equilibrium from
## scripts/equilibrium.R -- because the function is actually not
## terrifically smooth, it's not really the equilbrium seed rain, but
## close enough for our purposes.
p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- 505.55407
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast_control()) # A bit faster
p$set_control_parameters(list(schedule_verbose=TRUE))

schedule <- build_schedule(p)
schedule$use_ode_times <- TRUE

## Now, let's get a view of the fitness landscape with respect to one
## parameter.
lma <- p[[1]]$parameters[["lma"]]
n_mutants <- 11

## Add the mutants; they will all be introduced with a seed rain of 1,
## though they will not influence the environent at all.
## p_with_mutants <- p$copy()
## Add a vector of new lmas.  There will be n_mutants+1 mutant
## strategy, because we're adding the resident in as a mutant too for
## checking.
## p_with_mutants$add_strategy_mutant(new(Strategy, list(lma=lma)))
## for (i in seq(lma * 0.1, lma * 1.2, length.out=n_mutants)) {
##   p_with_mutants$add_strategy_mutant(new(Strategy, list(lma=i)))
## }
lma_mutant <- c(lma, seq(lma * 0.1, lma * 1.2, length.out=n_mutants))
p_with_mutants <- expand_parameters("lma", lma_mutant, p)
schedule_with_mutants <- expand_schedule(schedule, p_with_mutants$n_mutants)

## And run, producing only fitnesses.  The resident fitness should be
## unchanged here.
##
## Unfortunately this is very slow!
ebt_with_mutants <- run_ebt(p_with_mutants, schedule_with_mutants)
w_with_mutants <- ebt_with_mutants$fitnesses

lma_v <- sapply(seq_len(p_with_mutants$size),
                function(i) p_with_mutants[[i]]$parameters[["lma"]])

all.equal(lma_v, c(lma, lma_mutant))

##+ fitness_landscape
plot(lma_v[-(1:2)], w_with_mutants[-(1:2)], log="y")
points(lma_v[[1]], w_with_mutants[[1]] / p$seed_rain[[1]],
       col="red", pch=19)
points(lma_v[[2]], w_with_mutants[[2]],
       col="blue", pch=19, cex=.5)
