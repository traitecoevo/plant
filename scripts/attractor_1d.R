library(tree)
library(parallel)

## solve for 1D evolutionary attractor using root finding

p <- ebt_base_parameters()

root <- find_singularity_1D("lma", interval=c(0.075, 0.085),
                            p=p, tol=1e-04)

## check result by plotting selection gradient
xx <- seq(0.05, 0.1, length.out = 10)
yy <- selection_gradient("lma", xx, p)

plot(xx, yy, type = "b", xlab = "trait", ylab = "Selection gradient")
lines(spline(xx, yy, n=101), col="blue", lty=2)
abline(h = 0)
points(root$trait, 0, pch = 16)
