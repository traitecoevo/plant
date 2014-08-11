library(tree)

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$set_control_parameters(list(schedule_verbose=TRUE))

sys0 <- community(p0, "lma", seed_rain_initial=1e-3,
                  bounds=rbind(lma=c(0.01, 10)))

# Find point of maximum fitness
x <- tree:::max_fitness(sys0$trait_names, sys0$to_parameters(),
                        sys0$bounds)

# Make assembler
obj <- assembler_stochastic_naive(sys0, compute_viable_fitness=TRUE)

# Run assembler
set.seed(1)
obj$run_nsteps(7)

# Find 1D ESS
root <- find_singularity_1D(sys0$trait_names, interval = c(0.075, 0.085), p = sys0$to_parameters(), tol = 1e-04)

# Extract useful stuff
sys <- obj$get_community()
sys$to_parameters()
sys$to_schedule()
f <- sys$make_landscape()
m <- sys$traits(TRUE)
h <- obj$get_history()

# Extract stuff from history
extract_traits <- function(x){
	tree:::collect("traits", x, empty=NULL, loop=lapply, each=unlist, after=tree:::rbind_list)
	}

lapply(h,extract_traits)
