library(tree)

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$set_control_parameters(list(schedule_verbose=TRUE))

max_bounds <- rbind(lma=c(0.01, 10))

# Set parameters of ecological model

sys0 <- community(p0, "lma", seed_rain_initial=1e-3,
                  bounds=max_bounds)

# Find point of maximum fitness
# Todo: should this be a function of community?
# only need to do this is default value for trait has neg fitness
# Function has problems when bounds are too extreme. Use bracketing
# to isolate where peak is?
max_fit <- tree:::max_fitness(sys0$trait_names, sys0$to_parameters(),
                        sys0$bounds)



if(max_fit$objective < 0)
	stop("Finished: No region of positive fitness exists")


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
