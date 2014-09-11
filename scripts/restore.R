library(tree)

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast_control())
p0$set_control_parameters(list(schedule_verbose=TRUE))
p0$disturbance <- new(Disturbance, 10)

max_bounds <- rbind(lma=c(0.01, 10))

sys0 <- community(p0, "lma", bounds=max_bounds)
filename <- "restore.rds"
obj <- assembler_sample_positive(sys0, approximate_type="gp",
                                 compute_viable_fitness=TRUE,
                                 jump_to_attractor=FALSE,
                                 filename=filename)

set.seed(1)
obj$run_nsteps(2)

## This is something that needs dealing with, as it's a real pain.
obj$get_community()$add_approximate_landscape(type="gp")
obj$history[[length(obj$history)]] <- obj$get_community()$serialise()
## This is surprisingly slow.  That suggests we've still got a few R6
## objects causing grief here.
obj$save_to_file()

filename_save1 <- paste0(filename, ".1")
h1 <- add_approximate_landscapes(obj$history, p0,
                                 filename=filename_save1)

## Here are the landscapes over time:
h1 <- obj$history
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

## Super fast on the way back through, because history is already
## there.
h1 <- add_approximate_landscapes(h1, p0)


## This takes almost a second.
x <- last(h1)
saveRDS(x, "tmp1.rds")
## But this saves almost instantly
y <- x
y$landscape_approximate <- NULL
saveRDS(y, "tmp2.rds")

x <- last(h1)$landscape_approximate
saveRDS(x, "tmp1.rds")

## OK, that's too slow.  Let's do this again.
