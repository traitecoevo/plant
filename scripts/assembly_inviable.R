## With disturbances of 1, 2, 4, 8, see how the community starts up.
## Only the last is actually viable.
library(tree)

max_fitness <- tree:::max_fitness
fitness_landscape_approximate <- tree:::fitness_landscape_approximate
seq_log_range <- tree:::seq_log_range

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$set_control_parameters(list(schedule_verbose=TRUE))
p0$strategy_default <- new(Strategy, list(lma=10))

f <- function(disturbance, p0) {
  p <- p0$copy()
  p$disturbance <- new(Disturbance, disturbance)
  bounds <- rbind(lma=c(0.01, 10))
  sys0 <- community(p, "lma", seed_rain_initial=1e-3, bounds=bounds)
  ## Just for comparison:
  peak <- max_fitness(sys0$trait_names, p$copy(), bounds)
  ## Spin up the community:
  ok <- sys0$set_viable_bounds(find_max_if_negative=TRUE)
  fitness <- fitness_landscape_approximate(sys0)
  list(ok=ok, bounds=sys0$bounds, fitness=fitness, peak=peak)
}

res <- lapply(c(1, 2, 4, 8), f, p0)

par(mfrow=c(2, 2))
for (x in res) {
  lma <- seq_log_range(x$bounds, 200)
  plot(lma, x$fitness(lma), type="l", log="x")
  points(x$peak, attr(x$peak, "fitness"), pch=19)
}
