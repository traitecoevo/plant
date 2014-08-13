## Case getting the community load/save working.
library(tree)

ebt_base_parameters <- tree:::ebt_base_parameters
restore_community <- tree:::restore_community
fitness_landscape_approximate <- tree:::fitness_landscape_approximate
seq_log_range <- tree:::seq_log_range
restore_history <- tree:::restore_history

## This is from our successional_diversity project: generating new
## parameters from a few hyperparameters:
assembler_lma_parameters <- function(time_disturbance, slope) {
  ## Specialised parameters:
  model <- list(c_r1=0.5,
                c_r2=0,
                B4=slope,
                ## keeps mean centered at global mean:
                a4=10^(0.1369 + slope*-0.9819))
  ## Build up a parameters object from all of that
  p <- ebt_base_parameters()
  p$strategy_default <- new(Strategy, model)
  p$disturbance <- new(Disturbance, time_disturbance)
  p
}

## Start with a simple commiunity:
p0 <- assembler_lma_parameters(10, 1.7)
sys0 <- community(p0, "lma", bounds=rbind(lma=c(0.01, 10)))
## Add some species, and check that it serialises nicely:
sys0$add_traits(c(.1, 1, 10))
sys0$serialise()

## Dump tthis into a file and reload.
filename <- "test_restore.rds"
saveRDS(sys0$serialise(), filename)
obj <- restore_community(readRDS(filename), p0, FALSE)
## when computing the fitness landscape, ode times will be computed
## automatically.
f <- fitness_landscape_approximate(obj)

obj <- restore_community(dat, p0, TRUE)
f <- fitness_landscape_approximate(obj)

lma <- seq_log_range(obj$bounds, 400)
plot(lma, f(lma), log="x", type="l")

set.seed(1)
sys0 <- community(p0, "lma", bounds=rbind(lma=c(0.01, 10)))
obj <- assembler_stochastic_naive(sys0, compute_viable_fitness=TRUE)
obj$run_nsteps(5)

## Save the history to a file:
saveRDS(obj$get_history(), filename)
## Reload the history:
h0 <- restore_history(readRDS(filename), p0, FALSE)

x <- h0[[6]]
x$private$last_schedule # NULL
f <- fitness_landscape_approximate(x)
x$private$last_schedule # No longer NULL
h0[[6]] # also here

## Here's the fitness landscape:
lma <- seq_log_range(x$bounds, 400)
plot(lma, f(lma), log="x", type="l")
abline(v=x$traits(TRUE))
