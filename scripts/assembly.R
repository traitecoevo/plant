library(tree)

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$set_control_parameters(list(schedule_verbose=TRUE))

bounds_lma  <- c(0.01, 10)
bounds <- rbind(lma=bounds_lma)
colnames(bounds) <- c("lower", "upper")

sys0 <- community(p0, seed_rain_initial=1e-3)
obj <- assembler_stochastic_naive(sys0, bounds)
set.seed(1)

obj$nsteps(7)

sys <- obj$get_community()
sys$to_parameters()
sys$to_schedule()
f <- sys$make_landscape()

m <- sys$traits(TRUE)

h <- obj$get_history()

extract_traits <- function(x){
	tree:::collect("traits", x, empty=NULL, loop=lapply, each=unlist, after=tree:::rbind_list)
	}

lapply(h,extract_traits)
