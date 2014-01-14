source("helper-tree.R")

context("Disturbance")

## First, set things up the way that Daniel had it.
make.disturbance <- function(site.mean) {
  pow <- function(a, b) a^b
  psi <- 2.0
  ## solve lam as function of site average patch age
  lam <- pow(gamma(1.0/psi)/psi/site.mean, psi)
  ## solve for density age zero
  p0 <- psi*pow(lam, 1.0/psi)/gamma(1.0/psi)

  Pi <- function(age)
    exp(-lam*pow(age, psi))
  rate <- function(age)
    lam*psi*pow(age, psi-1)
  freq <- function(age)
    p0 * Pi(age)
  weight <- function(time.start, time)
    Pi(time)/ Pi(time.start)
  cdf <- function(p)
    pow(log(p) / -lam, 1/psi)

  list(site.mean=site.mean, psi=psi, lam=lam, p0=p0,
       rate=rate, freq=freq, Pi=Pi, weight=weight, cdf=cdf)
}

m <- 30.0
obj <- new(Disturbance, m)

test_that("Object has correct parameter", {
  expect_that(obj$mean_interval, is_identical_to(m))
})

disturbance <- make.disturbance(obj$mean_interval)

tt <- seq(0, 100, length=101)
test_that("Disturbance calculations are expected", {
  p.t <- sapply(tt, function(t) obj$pr_survival(t))
  expect_that(p.t, equals(disturbance$Pi(tt)))

  t.start <- 5
  p.t2 <- sapply(tt, function(t) obj$pr_survival_conditional(t, t.start))
  expect_that(p.t2, equals(disturbance$weight(t.start, tt)))

  ## Check of the conditional distribution approach:
  expect_that(p.t2,
              equals(p.t / obj$pr_survival(t.start)))

  expect_that(sapply(tt, function(t) obj$pr_survival_conditional(t, 0)),
              equals(p.t))
  expect_that(sapply(tt, function(t) obj$density(t)),
              equals(disturbance$freq(tt)))

  expect_that(disturbance$cdf(p.t), equals(tt))
  expect_that(sapply(p.t, function(p) obj$cdf(p)), equals(tt))
})

## Now, look at the rest of the issues.

## A Weibull distribution has PDF proportional to
##   (x/lambda)^(k-1) * exp(-(x/lambda)^k)
## None of the functions in the above set have that form, though the
## CDF is close.
scale <- with(disturbance, lam^(-1/psi))
shape <- disturbance$psi
expect_that(pweibull(tt, shape, scale, FALSE),
            equals(disturbance$weight(0, tt)))
