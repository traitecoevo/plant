source("helper-tree.R")

context("Gradient")

make_f <- function(n) {
  force(n)
  function(x) {
    if (is.matrix(x)) {
      if (ncol(x) != n) {
        stop("Wrong number of columns")
      }
    } else {
      if (length(x) != n) {
        stop("Wrong number of elements")
      }
      x <- rbind(x, deparse.level=0)
    }
    rowSums(rep(seq_len(n), each=nrow(x)) * (exp(x) - x)) / n
  }
}

make_g <- function(n) {
  vcv <- diag(seq_len(n))
  set.seed(1)
  vcv[upper.tri(vcv)] <- runif(n * (n - 1) / 2, -.5, .5)
  vcv[lower.tri(vcv)] <- t(vcv)[lower.tri(vcv)]
  function(x) {
    if (!(length(x) == n || ncol(x) == n)) {
      stop("Invalid size x")
    }
    mvtnorm::dmvnorm(x, sigma=vcv, log=TRUE)
  }
}

test_that("bootstrap", {
  sc2.f <- function(x){
    n <- length(x)
    sum((1:n) * (exp(x) - x)) / n
  }
  n <- 5
  f <- make_f(n)
  set.seed(1)
  x <- runif(n)
  expect_that(f(x), equals(sc2.f(x)))

  xx <- matrix(runif(n * 7), ncol=n)
  yy1 <- apply(xx, 1, sc2.f)
  yy2 <- apply(xx, 1, f)
  yy3 <- f(xx)
  expect_that(yy2, equals(yy1))
  expect_that(yy3, is_identical_to(yy2))
})

test_that("gradient (2d)", {
  n <- 2
  f <- make_f(n)
  set.seed(1)
  x <- runif(n)

  gr <- numDeriv::grad(f, x)
  expect_that(gradient(f, x), equals(gr))

  h <- numDeriv::hessian(f, x, method.args=list(d=0.0001))
  expect_that(Hessian(f, x), equals(h))

  obj <- slope_info(f, x)
  expect_that(obj$x,  equals(x))
  expect_that(obj$fx, equals(f(x)))
  expect_that(obj$gr, equals(gr))
  expect_that(obj$H,  equals(h))
})

test_that("gradient (2d quadratic)", {
  n <- 2
  f <- make_g(n)
  set.seed(1)
  x <- runif(n)

  gr <- numDeriv::grad(f, x)
  expect_that(gradient(f, x), equals(gr))

  h <- numDeriv::hessian(f, x, method.args=list(d=0.0001))
  expect_that(Hessian(f, x), equals(h))

  D <- numDeriv::genD(f, x)$D
  expect_that(D, equals(numDeriv::genD(f, x)$D))
  expect_that(D[seq_len(n)], equals(gr))
  expect_that(D[-seq_len(n)], equals(h[upper.tri(h, TRUE)]))

  obj <- slope_info(f, x)
  expect_that(obj$x,  equals(x))
  expect_that(obj$fx, equals(f(x)))
  expect_that(obj$gr, equals(gr))
  expect_that(obj$H,  equals(h))
})

test_that("gradient (5d)", {
  n <- 5
  f <- make_f(n)
  set.seed(1)
  x <- runif(n)

  gr <- numDeriv::grad(f, x)
  expect_that(gradient(f, x), equals(gr))

  h <- numDeriv::hessian(f, x, method.args=list(d=0.0001))
  expect_that(Hessian(f, x), equals(h))

  D <- bates_watts_D(f, x)$D
  expect_that(D, equals(numDeriv::genD(f, x)$D))
  expect_that(D[seq_len(n)], equals(gr))
  expect_that(D[-seq_len(n)], equals(h[upper.tri(h, TRUE)]))

  obj <- slope_info(f, x)
  expect_that(obj$x,  equals(x))
  expect_that(obj$fx, equals(f(x)))
  expect_that(obj$gr, equals(gr))
  expect_that(obj$H,  equals(h))
})

test_that("gradient (5D quadratic)", {
  n <- 5
  f <- make_g(n)
  set.seed(1)
  x <- runif(n)

  gr <- numDeriv::grad(f, x)
  expect_that(gradient(f, x), equals(gr))

  h <- numDeriv::hessian(f, x, method.args=list(d=0.0001))
  expect_that(Hessian(f, x), equals(h))

  D <- bates_watts_D(f, x)$D
  expect_that(D, equals(numDeriv::genD(f, x)$D))
  expect_that(D[seq_len(n)], equals(gr))
  expect_that(D[-seq_len(n)], equals(h[upper.tri(h, TRUE)]))

  obj <- slope_info(f, x)
  expect_that(obj$x,  equals(x))
  expect_that(obj$fx, equals(f(x)))
  expect_that(obj$gr, equals(gr))
  expect_that(obj$H,  equals(h))
})
