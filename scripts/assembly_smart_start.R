library(tree)

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$set_control_parameters(list(schedule_verbose=TRUE))
p0$disturbance <- new(Disturbance, 10)

max_bounds <- rbind(lma=c(0.01, 10))

# Set parameters of ecological model

## Get a community started, and take one step:
sys0 <- community(p0, "lma", bounds=max_bounds)
obj <- assembler_sample_positive(sys0)
obj$step()

## Do we have an approximate landscape for the empty community?
obj$history[[1]]

## Check the load/save cycle:
filename <- "test_restore.rds"
saveRDS(obj$history, filename)
tmp <- readRDS(filename)
f <- tmp[[1]]$landscape_approximate
lma <- tree:::seq_log_range(tmp[[1]]$bounds, 200)
plot(lma, f(lma), log="x", type="l")
h <- tree:::restore_history(tmp)

## Starting again...
sys0 <- community(p0, "lma", bounds=max_bounds)
sys0$set_viable_bounds()
sys0$jump_to_attractor()
f <- fitness_landscape_approximate(sys0)

lma <- seq_log_range(sys0$bounds, 400)
plot(lma, f(lma), type="l", log="x")
abline(v=sys0$traits(TRUE), h=0, col="red")

## And again, using an argument to the assembler:
sys0 <- community(p0, "lma", bounds=max_bounds)
obj <- assembler_sample_positive(sys0, compute_viable_fitness=TRUE,
                                 jump_to_attractor=TRUE)

vcv <- mutational_vcv_proportion(max_bounds, 0.001)
sys0 <- community(p0, "lma", bounds=max_bounds)
obj <- assembler_stochastic_naive(sys0, vcv,
                                  compute_viable_fitness=TRUE,
                                  jump_to_attractor=TRUE)
