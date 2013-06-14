source("helper-tree.R")

context("Disturbance")

## First, set things up the way that Daniel had it.
make.disturbance <- function(site.mean) {
  pow <- function(a, b) a^b
  psi <- 2.0
  ## solve lam as function of site average patch age
  lam <- pow(gamma(1.0/psi)/psi/site.mean, psi);
  ## solve for density age zero
  p0 <- psi*pow(lam, 1.0/psi)/gamma(1.0/psi);

  Pi <- function(age)
    exp(-lam*pow(age, psi))
  rate <- function(age)
    lam*psi*pow(age, psi-1)
  freq <- function(age)
    p0 * Pi(age)
  weight <- function(time.start, time)
    Pi(time)/ Pi(time.start)

  list(site.mean=site.mean, psi=psi, lam=lam, p0=p0,
       rate=rate, freq=freq, Pi=Pi, weight=weight)
}

obj <- new(Disturbance)

expected <- list(mean_disturbance_interval=30.0)
expect_that(obj$parameters, is_identical_to(expected))

new.p <- list(mean_disturbance_interval=20.0)
obj$set_parameters(new.p)
expect_that(obj$parameters,
            is_identical_to(modifyList(expected, new.p)))

disturbance <- make.disturbance(obj$parameters$mean_disturbance_interval)

tt <- seq(0, 100, length=101)
expect_that(sapply(tt, function(t) obj$survival_probability(0, t)),
            equals(disturbance$weight(0, tt)))

## Now, look at the rest of the issues.

## A Weibull distribution has PDF proportional to
##   (x/lambda)^(k-1) * exp(-(x/lambda)^k)
## None of the functions in the above set have that form, though the
## CDF is close.
scale <- with(disturbance, lam^(-1/psi))
shape <- disturbance$psi
expect_that(pweibull(tt, shape, scale, FALSE),
            equals(disturbance$weight(0, tt)))
