## This code is ported over from Rich's competition-kernels project.

library(tree)
library(parallel)

p <- ebt_base_parameters()

## Find the viable range for lma for the default parameters:
trait <- "lma"
r <- tree:::viable_fitness(trait, p, bounds=c(1E-5, 1E3), log_scale=TRUE)

## Construct a fitness landscape across this region showing it is the
## right set of bounds:
values <- seq_log(r[[1]], r[[2]], length=51)
w <- tree:::max_growth_rate(trait, values, p)
plot(w ~ values, log="x", type="l", ylim=range(0, w),
     xlab=trait, ylab="Max growth rate")
abline(h=0, lty=2, col="red")
