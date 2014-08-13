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
obj$get_history()[[1]]

## Check the load/save cycle:
filename <- "test_restore.rds"
saveRDS(obj$get_history(), filename)
tmp <- readRDS(filename)
f <- tmp[[1]]$landscape_approximate
lma <- tree:::seq_log_range(tmp[[1]]$bounds, 200)
plot(lma, f(lma), log="x", type="l")
foo <- tree:::restore_history(tmp, p0)
foo$landscape_approximate

## Starting again...
sys0 <- community(p0, "lma", bounds=max_bounds)
obj <- assembler_sample_positive(sys0, compute_viable_fitness=TRUE)

# solve for 1D ESS
# todo: should this be a function within assembler? Could then save result directly in community, icnluding seed rain
root <- find_singularity_1D(sys0$trait_names, interval = obj$get_community()$bounds, p = sys0$to_parameters(), tol = 1e-03)
root <- list(root= 0.07989016, f.root=-0.02790277, iter=16, estim.prec=5e-05)

#run landscape and check for positive areas

# Run assembler using above solution as starting point

set.seed(1)
sys1 <- community(p0, "lma", seed_rain_initial=1e-3, bounds=obj$get_community()$bounds)
sys1$add_traits(rbind(lma=root$root))
obj2 <- assembler_sample_positive(sys1, bounds <- obj$get_community()$bounds, compute_viable_fitness=FALSE)
obj2$step("to_equilibrium") # This not working
