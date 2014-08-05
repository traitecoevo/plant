library(tree)
library(parallel)

## solve for 1D evolutionary attractor using root finding

p <- new(Parameters)
p$set_parameters(list(patch_area = 1))  # See issue #13
p$set_control_parameters(fast_control())  # A bit faster

root <- find_singularity_1D("lma", interval = c(0.075, 0.085), p = p, tol = 1e-04)

# check result by plotting selection gradient

xx <- seq(0.05, 0.1, length.out = 10)
yy <- selection_gradient("lma", xx, p)

plot(xx, yy, type = "b", xlab = "trait", ylab = "Selection gradient")
abline(h = 0)
points(root$root, root$f.root, pch = 16)
