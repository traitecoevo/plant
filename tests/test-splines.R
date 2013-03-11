source("helper-tree.R")

context("Splines")

target <- function(x) sin(2*x)

xx <- seq(0, 2*pi, length=101)
yy <- target(xx)

## Basic splines:
s <- new(Spline)
s$init(xx, yy)
expect_that(s$xy(),
            is_identical_to(cbind(xx, yy, deparse.level=0)))

## Test data:
xx.cmp <- seq(min(xx), max(xx), length=2*length(xx))
yy.cmp <- target(xx.cmp)

yy.R <- spline(xx, yy, xout=xx.cmp)$y
yy.C <- s$eval(xx.cmp)

expect_that(yy.R, equals(yy.cmp, tolerance=1e-6))
expect_that(yy.C, equals(yy.cmp, tolerance=1e-6))

## Adaptive Splines:

a <- new(AdaptiveSplineR, target, new.env(), min(xx), max(xx))
expect_that(nrow(a$xy()), equals(241))

xx.eval <- a$xy()[,1]
expect_that(a$xy()[,2], is_identical_to(target(xx.eval)))

xx.mid <- (xx.eval[-1] + xx.eval[-length(xx.eval)]) / 2
yy.mid <- target(xx.mid)
zz.mid <- a$eval(xx.mid)
expect_that(zz.mid, equals(yy.mid, tolerance=2e-8))

## TODO: Use tolerance from object:
err <- pmax(abs(zz.mid - yy.mid), abs(1 - zz.mid / yy.mid))
expect_that(all(err < 1e-6), is_true())
