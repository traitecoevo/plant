library(tree)
library(siefecor)

max_fitness <- tree:::max_fitness
mutational_vcv_proportion <- tree:::mutational_vcv_proportion
collect <- tree:::collect
rbind_list <- tree:::rbind_list

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$set_control_parameters(list(schedule_verbose=TRUE))
p0$strategy_default <- new(Strategy, list(lma=10))
p0$disturbance <- new(Disturbance, 7)

max_bounds <- rbind(lma=c(0.01, 10))

sys0 <- community(p0, "lma", seed_rain_initial=1e-3,
                  bounds=max_bounds)
ok <- sys0$set_viable_bounds(find_max_if_negative=TRUE)

f <- fitness_landscape_approximate(sys0)
lma <- seq_log_range(sys0$bounds, 400)
plot(f(lma) ~ lma, type="l", log="x")

g <- fitness_landscape_approximate(sys0, "gp")
plot(f(lma) ~ lma, type="l", log="x")
lines(g(lma) ~ lma, type="l", col="red")

sys0 <- community(p0, "lma", bounds=max_bounds)
obj <- assembler_sample_positive(sys0, approximate_type="gp")






## OK, so the approximate landscape code is the next target.

## Currently we do fitness_landscape_approximate, which takes a few
## points for us and returns a spline function.
##
## This is used within `make_births_sample_positive`.
##
## The simplest thing seems to be to set up an argument-free version
## and use that.  We don't really have support for a general version
## of the non-gp approach yet either.


