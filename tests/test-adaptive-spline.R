source("helper-tree.R")

context("AdaptiveSpline")

target <- function(x) sin(2*x)
r <- c(0, 2*pi)
a <- tree_module$test_adaptive_spline(target, new.env(), r[1], r[2])

expect_that(nrow(a$xy), equals(241))

xx.eval <- a$xy[,1]
expect_that(a$xy[,2], is_identical_to(target(xx.eval)))

xx.mid <- (xx.eval[-1] + xx.eval[-length(xx.eval)]) / 2
yy.mid <- target(xx.mid)
zz.mid <- a$eval(xx.mid)
expect_that(zz.mid, equals(yy.mid, tolerance=2e-8))

err <- pmax(abs(zz.mid - yy.mid), abs(1 - zz.mid / yy.mid))
expect_that(all(err < 1e-6), is_true())
