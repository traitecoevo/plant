source("helper-tree.R")

context("Spline")

target <- function(x) sin(2*x)

xx <- seq(0, 2*pi, length=101)
yy <- target(xx)

## Basic splines:
s <- new(Spline)
expect_that(s$min, equals(Inf))
expect_that(s$max, equals(-Inf))

s$init(xx, yy)
expect_that(s$xy,
            is_identical_to(cbind(xx, yy, deparse.level=0)))
expect_that(c(s$min, s$max), is_identical_to(range(xx)))

## Test data:
xx.cmp <- seq(min(xx), max(xx), length=2*length(xx))
yy.cmp <- target(xx.cmp)

yy.R <- spline(xx, yy, xout=xx.cmp)$y
yy.C <- s$eval(xx.cmp)

expect_that(yy.R, equals(yy.cmp, tolerance=1e-6))
expect_that(yy.C, equals(yy.cmp, tolerance=1e-6))
