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

  list(site.mean=site.mean, psi=psi, lam=lam, p0=p0,
       rate=rate, freq=freq, Pi=Pi, weight=weight)
}

test_that("Creation", {
  m <- 30.0
  obj <- Disturbance(m)
  expect_is(obj, "Disturbance")
  expect_is(obj, "R6")
  expect_identical(obj$mean_interval, m)
})

test_that("Disturbance calculations are expected", {
  m <- 30.0
  obj <- Disturbance(m)

  disturbance <- make_disturbance(obj$mean_interval)
  tt <- seq(0, 100, length.out=101)

  p_t <- sapply(tt, function(t) obj$pr_survival(t))
  expect_equal(p_t, disturbance$Pi(tt))

  t_start <- 5
  p_t2 <- sapply(tt, function(t) obj$pr_survival_conditional(t, t_start))
  expect_equal(p_t2, disturbance$weight(t_start, tt))

  ## Check of the conditional distribution approach:
  expect_equal(p_t2, p_t / obj$pr_survival(t_start))

  expect_equal(sapply(tt, function(t) obj$pr_survival_conditional(t, 0)), p_t)

  ## density is vectorised
  expect_equal(obj$density(tt), disturbance$freq(tt))

  ## Now, look at the rest of the issues.

  ## A Weibull distribution has PDF proportional to
  ##   (x/lambda)^(k-1) * exp(-(x/lambda)^k)
  ## None of the functions in the above set have that form, though the
  ## CDF is close.
  scale <- with(disturbance, lam^(-1/psi))
  shape <- disturbance$psi
  expect_equal(pweibull(tt, shape, scale, FALSE), disturbance$weight(0, tt))
})

# test_that("Reference survival eps gives correct running time", {
#   ## From falster-traitdiversity: src/base/ebt/site.cpp,
#   ## site::solve_patchage_dist():
#   f1 <- function(mean) {
#     2.633 * mean / 3.0 * 4.0
#   }
#
#   ## Wrapper for getting same out of our disturbance class:
#   ctrl <- Control()
#   eps <- ctrl$schedule_patch_survival
#   f2 <- function(mean) {
#     Disturbance(mean)$cdf(eps)
#   }
#
#   age <- seq(1, 200, length.out=101)
#   y1 <- f1(age)
#   y2 <- sapply(age, f2)
#
#   expect_equal(y1, y2)
# })
