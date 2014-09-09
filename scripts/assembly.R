library(tree)

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

# Find point of maximum fitness
x <- max_fitness(sys0$trait_names, sys0$to_parameters(), sys0$bounds)

# Make assembler
vcv <- mutational_vcv_proportion(max_bounds, 0.001)

sys0 <- community(p0, "lma", seed_rain_initial=1e-3,
                  bounds=max_bounds)

obj_n <- assembler_sample_positive(sys0, "naive")
obj_g <- assembler_sample_positive(sys0, "gp")

obj <- assembler_stochastic_naive(sys0, vcv)

# Run assembler
set.seed(1)
obj$run_nsteps(7)
