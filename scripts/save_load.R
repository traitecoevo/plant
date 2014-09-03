## Case getting the community load/save working.
library(tree)

ebt_base_parameters <- tree:::ebt_base_parameters
restore_community <- tree:::restore_community
fitness_landscape_approximate <- tree:::fitness_landscape_approximate
seq_log_range <- tree:::seq_log_range
restore_history <- tree:::restore_history
mutational_vcv_proportion <- tree:::mutational_vcv_proportion

## This is from our successional_diversity project: generating new
## parameters from a few hyperparameters:
assembler_lma_parameters <- function(time_disturbance, slope) {
  ## Specialised parameters:
  model <- list(c_r1=0.5,
                c_r2=0,
                B4=slope)
  ## Build up a parameters object from all of that
  p <- ebt_base_parameters()
  p$strategy_default <- new(Strategy, model)
  p$disturbance <- new(Disturbance, time_disturbance)
  p
}

max_bounds <- rbind(lma=c(0.01, 10))

## Start with a simple commiunity:
p0 <- assembler_lma_parameters(10, 1.7)
sys0 <- community(p0, "lma", bounds=max_bounds)
## Add some species, and check that it serialises nicely:
sys0$add_traits(c(.1, 1, 10))
sys0$serialise()

## Dump this into a file and reload.
filename <- "test_restore.rds"
saveRDS(sys0$serialise(), filename)
obj <- restore_community(readRDS(filename))
## when computing the fitness landscape, ode times will be computed
## automatically.
f <- fitness_landscape_approximate(obj)

lma <- seq_log_range(obj$bounds, 400)
plot(lma, f(lma), log="x", type="l")

vcv <- mutational_vcv_proportion(max_bounds, 0.001)
set.seed(1)
sys0 <- community(p0, "lma", bounds=rbind(lma=c(0.01, 10)))
obj <- assembler_stochastic_naive(sys0, vcv)
obj$run_nsteps(5)

## Save the history to a file:
saveRDS(obj$get_history(), filename)
## Reload the history:
h0 <- restore_history(readRDS(filename), recompute=FALSE)

x <- h0[[6]]
x$private$last_schedule # NULL
f <- fitness_landscape_approximate(x)
x$private$last_schedule # No longer NULL
h0[[6]]$private$last_schedule # also here

## Here's the fitness landscape:
lma <- seq_log_range(x$bounds, 400)
plot(lma, f(lma), log="x", type="l")
abline(v=x$traits(TRUE))

## Then again with the sampling version:
set.seed(1)
sys0 <- community(p0, "lma", bounds=rbind(lma=c(0.01, 10)))
obj <- assembler_sample_positive(sys0)
obj$run_nsteps(5, "to_equilibrium")
