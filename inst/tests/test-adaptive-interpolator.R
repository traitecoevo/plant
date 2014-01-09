source("helper-tree.R")

context("AdaptiveInterpolator")

target <- function(x) sin(2*x)
r <- c(0, 2*pi)
s <- test_adaptive_interpolator(target, new.env(), r[1], r[2], FALSE)

test_that("Interpolator is the size expected", {
  expect_that(s$size, equals(241))
  expect_that(nrow(s$xy), equals(s$size))
})

xx.eval <- s$xy[,1]
test_that("Evaluated points are exactly known", {
  expect_that(s$xy[,2], is_identical_to(target(xx.eval)))
})

test_that("Tolerance is as expected", {
  xx.mid <- (xx.eval[-1] + xx.eval[-length(xx.eval)]) / 2
  yy.mid <- target(xx.mid)
  zz.mid <- s$eval(xx.mid)
  expect_that(zz.mid, equals(yy.mid, tolerance=2e-8))
  err <- pmax(abs(zz.mid - yy.mid), abs(1 - zz.mid / yy.mid))
  expect_that(all(err < 1e-6), is_true())
})

## And again with Akima splines:

source("helper-tree.R")

a <- test_adaptive_interpolator(target, new.env(), r[1], r[2], TRUE)

test_that("Interpolator is the size expected", {
  expect_that(a$size, equals(505))
  expect_that(nrow(a$xy), equals(a$size))
})

xx.eval <- a$xy[,1]
test_that("Evaluated points are exactly known", {
  expect_that(a$xy[,2], is_identical_to(target(xx.eval)))
})

test_that("Tolerance is as expected", {
  xx.mid <- (xx.eval[-1] + xx.eval[-length(xx.eval)]) / 2
  yy.mid <- target(xx.mid)
  zz.mid <- a$eval(xx.mid)
  ## Looks like we bottomed out in terms of refinement:
  ## Akima splines aren't going to be much use.
  expect_that(zz.mid, equals(yy.mid, tolerance=6e-7))
  err <- pmax(abs(zz.mid - yy.mid), abs(1 - zz.mid / yy.mid))
  expect_that(all(err < 4e-5), is_true())
})
