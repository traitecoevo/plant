source("helper-tree.R")

context("Fake light environment")

## Need a little 3d function; doesn't need to be crazy; this is going
## to be BM with drift.

target <- function(t, x)
  1 - dnorm(x, t * 0.01, 0.1 + t/10)

## Times:
tt <- seq(0, 5, length.out=11)

## Start with the same x points (then do overlapping sets then do
## extrapolation)
interpolator <- function(x, y) {
  s <- new(Interpolator)
  s$init(x, y)
  s
}

## Easiest case: the splines have the same x (height) spacing so
## extrapolation and extra point calculation is not an issue.
test_that("Shared x spacing", {
  xx <- seq(-5, 5, length.out=51)
  ss <- lapply(tt, function(t) interpolator(xx, target(t, xx)))
  obj <- new(FakeLightEnvironment, tt, ss)
  expect_that(obj$current$xy, is_identical_to(ss[[1]]$xy))

  ## Setting to the first, second and last times all return the
  ## equilvant splines exactly.
  obj$set_time(tt[[1]])
  expect_that(obj$current$xy, is_identical_to(ss[[1]]$xy))

  obj$set_time(tt[[2]])
  expect_that(obj$current$xy, is_identical_to(ss[[2]]$xy))

  i <- length(tt)
  obj$set_time(tt[[i]])
  expect_that(obj$current$xy, is_identical_to(ss[[i]]$xy))

  ## Then try merging splines in the middle:
  i <- 5:6
  i1 <- i[1]
  i2 <- i[2]
  obj$set_time(mean(tt[i]))
  expect_that(obj$current$xy,
              equals((ss[[i1]]$xy + ss[[i2]]$xy)/2))

  ## and then at some random point
  set.seed(1)
  t <- runif(1, tt[i1], tt[i2])
  obj$set_time(t)

  w <- (t - tt[i2]) / (tt[i1] - tt[i2])
  z.approx <- ss[[i1]]$y * w + ss[[i2]]$y * (1-w)
  expect_that(obj$current$y, equals(z.approx))

  ## Some simple bounds checking:
  expect_that(obj$set_time(-1),          throws_error())
  expect_that(obj$set_time(max(tt) + 1), throws_error())
})

test_that("Different x spacings, shared end points", {
  r <- c(-5, 5)
  np <- 49
  set.seed(1)
  xx <- lapply(seq_along(tt), function(i)
               sort(c(r, runif(rpois(1, 49), r[1], r[2]))))
  ss <- lapply(seq_along(tt), function(i)
               interpolator(xx[[i]], target(tt[i], xx[[i]])))

  obj <- new(FakeLightEnvironment, tt, ss)
  expect_that(obj$current$xy, is_identical_to(ss[[1]]$xy))

  i <- 5:6
  i1 <- i[1]
  i2 <- i[2]
  set.seed(1)
  t <- runif(1, tt[i1], tt[i2])
  obj$set_time(t)

  w <- (t - tt[i2]) / (tt[i1] - tt[i2])
  xx.cmp <- sort(unique(c(ss[[i1]]$x, ss[[i2]]$x)))
  y1.cmp <- ss[[i1]]$eval(xx.cmp)
  y2.cmp <- ss[[i2]]$eval(xx.cmp)
  yy.cmp <- y1.cmp * w + y2.cmp * (1-w)

  expect_that(obj$current$size, equals(length(xx.cmp)))
  expect_that(obj$current$x,    equals(xx.cmp))
  expect_that(obj$current$y,    equals(yy.cmp))
})

test_that("Different x spacings, different end points", {
  r <- c(-5, 5)
  np <- 49
  set.seed(1)
  xx <- lapply(seq_along(tt), function(i)
               sort(c(r[1], runif(rpois(1, 49), r[1], r[2]))))
  xx <- xx[order(sapply(xx, max))]
  ss <- lapply(seq_along(tt), function(i)
               interpolator(xx[[i]], target(tt[i], xx[[i]])))

  obj <- new(FakeLightEnvironment, tt, ss)
  expect_that(obj$current$xy, is_identical_to(ss[[1]]$xy))

  i <- 5:6
  i1 <- i[1]
  i2 <- i[2]
  set.seed(1)
  t <- runif(1, tt[i1], tt[i2])
  obj$set_time(t)

  w <- (t - tt[i2]) / (tt[i1] - tt[i2])
  xx.cmp <- sort(unique(c(ss[[i1]]$x, ss[[i2]]$x)))
  y1.cmp <- c(ss[[i1]]$eval(xx.cmp[xx.cmp <= max(xx[[1]])]),
              rep(1, sum(xx.cmp > max(xx[[1]]))))
  y2.cmp <- ss[[i2]]$eval(xx.cmp)
  yy.cmp <- y1.cmp * w + y2.cmp * (1-w)

  expect_that(obj$current$size, equals(length(xx.cmp)))
  expect_that(obj$current$x,    equals(xx.cmp))
  expect_that(obj$current$y,    equals(yy.cmp))
})

test_that("Can't get spline past end time", {

})
