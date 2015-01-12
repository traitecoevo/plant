if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("Disturbance")

## First, set things up the way that Daniel had it.
make_disturbance <- function(site.mean) {
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

test_that("Creation", {
  m <- 30.0
  obj <- Disturbance(m)
  expect_that(obj, is_a("Disturbance"))
  expect_that(obj, is_a("R6"))
  expect_that(obj$mean_interval, is_identical_to(m))
})

test_that("Disturbance calculations are expected", {
  m <- 30.0
  obj <- Disturbance(m)

  disturbance <- make_disturbance(obj$mean_interval)
  tt <- seq(0, 100, length=101)

  p_t <- sapply(tt, function(t) obj$pr_survival(t))
  expect_that(p_t, equals(disturbance$Pi(tt)))

  t_start <- 5
  p_t2 <- sapply(tt, function(t) obj$pr_survival_conditional(t, t_start))
  expect_that(p_t2, equals(disturbance$weight(t_start, tt)))

  ## Check of the conditional distribution approach:
  expect_that(p_t2,
              equals(p_t / obj$pr_survival(t_start)))

  expect_that(sapply(tt, function(t) obj$pr_survival_conditional(t, 0)),
              equals(p_t))
  expect_that(sapply(tt, function(t) obj$density(t)),
              equals(disturbance$freq(tt)))

  expect_that(disturbance$cdf(p_t), equals(tt))
  expect_that(sapply(p_t, function(p) obj$cdf(p)), equals(tt))

  ## Now, look at the rest of the issues.

  ## A Weibull distribution has PDF proportional to
  ##   (x/lambda)^(k-1) * exp(-(x/lambda)^k)
  ## None of the functions in the above set have that form, though the
  ## CDF is close.
  scale <- with(disturbance, lam^(-1/psi))
  shape <- disturbance$psi
  expect_that(pweibull(tt, shape, scale, FALSE),
              equals(disturbance$weight(0, tt)))
})

test_that("Reference survival eps gives correct running time", {
  ## From falster-traitdiversity: src/base/ebt/site.cpp,
  ## site::solve_patchage_dist():
  f1 <- function(mean) {
    2.633 * mean / 3.0 * 4.0
  }

  ## Wrapper for getting same out of our disturbance class:
  ctrl <- Control()
  eps <- ctrl$schedule_default_patch_survival
  f2 <- function(mean) {
    Disturbance(mean)$cdf(eps)
  }

  age <- seq(1, 200, length.out=101)
  y1 <- f1(age)
  y2 <- sapply(age, f2)

  expect_that(y1, equals(y2))
})
