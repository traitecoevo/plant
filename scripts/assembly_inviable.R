## With disturbances of 1, 2, 4, 8, see how the community starts up.
## Only the last is actually viable.
library(tree)

max_fitness <- tree:::max_fitness
fitness_landscape_approximate <- tree:::fitness_landscape_approximate
mutational_vcv_proportion <- tree:::mutational_vcv_proportion
seq_log_range <- tree:::seq_log_range

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$set_control_parameters(list(schedule_verbose=TRUE))
p0$strategy_default <- new(Strategy, list(lma=10))

max_bounds <- rbind(lma=c(0.01, 10))

f <- function(disturbance, p0) {
  p <- p0$copy()
  p$disturbance <- new(Disturbance, disturbance)
  sys0 <- community(p, "lma", bounds=max_bounds)
  ## Just for comparison:
  peak <- max_fitness(sys0$trait_names, p$copy(), max_bounds)
  ## Spin up the community:
  ok <- sys0$set_viable_bounds(find_max_if_negative=TRUE)
  fitness <- fitness_landscape_approximate(sys0, bounds=max_bounds)
  list(ok=ok, bounds=sys0$bounds, fitness=fitness, peak=peak)
}

res <- lapply(c(1, 2, 4, 8), f, p0)

par(mfrow=c(2, 2))
for (x in res) {
  lma <- seq_log_range(if (is.null(x$bounds)) max_bounds else x$bounds, 200)
  plot(lma, x$fitness(lma), type="l", log="x")
  points(x$peak, attr(x$peak, "fitness"), pch=19)
}

## Should be nice and quick with an inviable case:
p <- p0$copy()
p$disturbance <- new(Disturbance, 2)
sys0 <- community(p, "lma", bounds=max_bounds)

vcv <- mutational_vcv_proportion(max_bounds, 0.001)
obj_sample <- assembler_sample_positive(sys0)
obj_sample$run_nsteps(20)
obj_sample$run_nsteps(20, "to_equilibrium")

obj_naive <- assembler_stochastic_naive(sys0, vcv)
obj_naive$run_nsteps(20)
obj_naive$run_nsteps(20, "to_equilibrium")
