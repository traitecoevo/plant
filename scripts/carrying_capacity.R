## This code is ported over from Rich's competition-kernels project.

library(tree)
library(parallel)

p <- new(Parameters)
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast_control()) # A bit faster
p$set_control_parameters(equilibrium_verbose())

## Find the viable range for lma for the default parameters:
trait <- "lma"
r <- tree:::viable_fitness(trait, p)

## Construct a fitness landscape across this region showing it is the
## right set of bounds:
values <- seq_log(r[[1]], r[[2]], length=51)
w <- tree:::max_growth_rate(trait, values, p)
plot(w ~ values, log="x", type="l", ylim=range(0, w),
     xlab=trait, ylab="Max growth rate")
abline(h=0, lty=2, col="red")
