source("helper-tree.R")

context("MultiSpline")

target <- function(x) c(sin(2*x), sin(3*x), sin(4*x))

xx <- seq(0, 2*pi, length=101)
yy <- t(sapply(xx, target))

ss <- new(MultiSpline, 3)
expect_that(ss$dim, equals(3))

ss$init(xx, yy)
expect_that(ss$size, equals(length(xx)))

## Test data:
xx.cmp <- seq(min(xx), max(xx), length=2*length(xx))
yy.cmp <- t(sapply(xx.cmp, target))

yy.R <- apply(yy, 2, function(yi)
              spline(xx, yi, xout=xx.cmp)$y)
yy.C <- ss$eval(xx.cmp)

expect_that(yy.R, equals(yy.cmp, tolerance=3e-6))
expect_that(yy.C, equals(yy.cmp, tolerance=3e-6))

expect_that(ss$x, equals(xx))
expect_that(ss$y, equals(yy))

ss$reset()

for ( i in seq_along(xx) )
  ss$add_point(xx[[i]], yy[i,])
ss$init_self()
expect_that(ss$size, equals(length(xx)))

expect_that(ss$eval(xx.cmp), is_identical_to(yy.C))
expect_that(ss$x, equals(xx))
expect_that(ss$y, equals(yy))
