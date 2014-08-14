library(tree)

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$set_control_parameters(list(schedule_verbose=TRUE))
p0$disturbance <- new(Disturbance, 10)

max_bounds <- rbind(lma=c(0.01, 10))

vcv <- mutational_vcv_proportion(max_bounds, 0.001)
sys0 <- community(p0, "lma", bounds=max_bounds)
filename <- "restore.rds"
obj <- assembler_stochastic_naive(sys0, vcv,
                                  compute_viable_fitness=TRUE,
                                  jump_to_attractor=FALSE,
                                  filename=filename)
set.seed(1)
obj$run_nsteps(5)

filename_save1 <- paste0(filename, ".1")
h1 <- add_approximate_landscapes(obj$get_history(), p0, filename_save1)

## Here are the landscapes over time:
lma <- seq_log_range(h1[[1]]$bounds, 400)
w <- sapply(h1, function(x) x$landscape_approximate(lma))
matplot(lma, w, type="l", lty=1, col="black", log="x")

## Then try with the saved history:
dat <- readRDS(filename)
filename_save2 <- paste0(filename, ".2")
h2 <- add_approximate_landscapes(dat, p0, filename_save2)

w2 <- sapply(h2, function(x) x$landscape_approximate(lma))
all.equal(w2, w)

## Try reloading both from files:
tmp1 <- readRDS("restore.rds.1")
tmp2 <- readRDS("restore.rds.2")
all.equal(sapply(tmp1, function(x) x$landscape_approximate(lma)), w)
all.equal(sapply(tmp2, function(x) x$landscape_approximate(lma)), w)
