library(tree)
source("assembly-analysis-fun.R")
source("assembly_lma-fun.R")
p <- assembler_lma_parameters(20, 1.7)
dat <- readRDS("from_daniel/3cc5589/naive-single-FALSE-20.0000-1.7000.rds")
dat <- readRDS("from_daniel/3cc5589/naive-to_equilibrium-TRUE-20.0000-1.7000.rds")

filename <- "for_presentation.rds"
if (!file.exists(filename)) {
  com <- restore_community(dat[[length(dat)]], p)
  add_approximate_landscape(com)
  saveRDS(com$serialise(), filename)
} else {
  com <- restore_community(readRDS(filename), p)
}

## TODO: Need to understand why adding a few points here is doing so
## appallingly.  Getting very different fitness values out.  This is
## in computing the second rounr of points

f <- com$make_landscape()
fa <- com$landscape_approximate

lma <- seq_log_range(com$bounds, 400)
plot(lma, fa(lma), type="l", log="x")

xy <- with(environment(com$landscape_approximate), data.frame(x, y))
points(xy)

xx <- seq_log_range(range(xy$x[xy$y > -1.5]), 50)
yy <- f(xx)

xy2 <- rbind(xy, data.frame(x=xx, y=yy))
xy2 <- xy2[order(xy2$x),]

fa2 <- tree:::splinefun_log(xy2$x, xy2$y)

plot(lma, fa2(lma), type="l", log="x")
lines(lma, fa(lma), col="grey")

## OK, that's not working nicely, try this:
com <- restore_community(dat[[length(dat)]], p)
xy3 <- tree:::fitness_landscape_grid(com, 150, bounds=com$bounds)

## And boom: second run through is broken.  
w <- com$make_landscape()
z <- w(xx[[1]])

x_use <- seq_log_range(range(xx), nrow(xy3) * 4 - 1)
y_use <- tree:::splinefun_log(xy3[,1], xy3[,2])(x_use)

plot(x_use, y_use, log="x", type="l")
points(xy3)

write.table(cbind(x_use, y_use), "for_presentation.txt",
            row.names=FALSE, col.names=FALSE)


######################################################################

library(tree)
source("assembly-analysis-fun.R")
source("assembly_lma-fun.R")
p <- assembler_lma_parameters(20, 1.7)
dat <- readRDS("from_daniel/3cc5589/naive-single-FALSE-20.0000-1.7000.rds")
dat <- readRDS("from_daniel/3cc5589/naive-to_equilibrium-TRUE-20.0000-1.7000.rds")

com <- restore_community(dat[[length(dat)]], p)
xy <- tree:::fitness_landscape_grid(com, 400, bounds=com$bounds)
write.table(xy, "for_presentation.txt",
            row.names=FALSE, col.names=FALSE)

######################################################################

library(siefecor)
filename <- "for_presentation.txt"
filename_log <- "for_presentation_log.txt"

dat <- siefecor:::load_locs_and_vals(filename)
dat[["locs"]] <- log(dat[["locs"]])
siefecor:::save_locs_and_vals(dat$locs, dat$vals, filename_log)

fit     <- siefecor::interp1d_files_helper(filename, frac_end=.15)
fit_log <- siefecor::interp1d_files_helper(filename_log, frac_end=.15)

plot(vals ~ locs, dat, type="l")
lines(vals ~ log(locs), fit, col="red")
lines(vals ~ locs, fit_log, col="blue")

fit     <- siefecor::interp1d_files_helper(filename,
                                           frac_start=.03, frac_end=.08)
fit_log <- siefecor::interp1d_files_helper(filename_log,
                                           frac_start=.03, frac_end=.08)

pdf("gp_50points.pdf")
plot(vals ~ locs, dat, type="l")
lines(vals ~ log(locs), fit, col="red")
lines(vals ~ locs, fit_log, col="blue")
dev.off()

## Generate a bunch of randoms:
res <- replicate(20,
                 siefecor::interp1d_files_helper(filename_log,
                                                 frac_start=.03,
                                                 frac_end=.08),
                 simplify=FALSE)

pdf("gp_variation.pdf")
plot(vals ~ locs, dat, type="l")
for (i in res) lines(vals ~ locs, i, col="#ff000055")
dev.off()
