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

sys <- obj$get_community()

foo <- tree:::landscape_approximate(sys, "gp")
x <- seq_log_range(sys$bounds, 50)
plot(x, foo$predict(x), log="x", type="l")

saveRDS(foo, file="tmp.rds")
bar <- readRDS("tmp.rds")

bar$predict(x)
bar$restore()
bar$predict(x)
