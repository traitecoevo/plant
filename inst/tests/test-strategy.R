source("helper-tree.R")

context("Strategy")

s <- new(Strategy)

obj <- s$parameters

expected <- list(
  B1     = 0.306,
  B4     = 1.71,
  B5     = 0,
  B6     = 0,
  Y      = 0.7,
  a1     = 5.44,
  a3     = 0.07,
  b      = 0.17,
  c_Rb   = 8024,
  c_Rl   = 21000,
  c_Rr   = 217,
  c_Rs   = 4012,
  c_acc  = 4,
  c_bio  = 0.0245,
  c_d0   = 0.01,
  c_d1   = 0.0,
  c_d2   = 5.5,
  c_d3   = 20,
  c_p1   = 150.36,
  c_p2   = 0.19,
  c_r1   = 1,
  c_r2   = 50,
  c_s0   = 0.1,
  eta    = 12,
  hmat   = 16.5958691,
  k_b    = 0.2,
  k_l0   = 0.4565855,
  k_r    = 1,
  k_s0   = 0.2,
  lma    = 0.1978791,
  n_area = 0.00187,
  rho    = 608,
  s      = 3.8e-5,
  theta  = 4669)

keys <- sort(names(expected))
expect_that(sort(names(obj)),
            is_identical_to(sort(names(expected))))
expect_that(obj[keys], equals(expected[keys]))

## Add some new parameters:
new1 <- list(c_acc=4.1, lma=1.2)
s$set_parameters(new1)
expect_that(s$parameters[names(new1)], is_identical_to(new1))
expect_that(s$parameters, equals(modifyList(expected, new1)))

## Generate a failure:
new2 <- list(unknown_key=1)
expect_that(s$set_parameters(new2), throws_error())

## And have the list remain unchanged and valid
expect_that(s$parameters, equals(modifyList(expected, new1)))

## Check that even if some elements are unknown, known parameters are
## not changed.
new3 <- list(Y=0.9, unknown_key=1, hmat=12)
expect_that(s$set_parameters(new3), throws_error())
expect_that(s$parameters, equals(modifyList(expected, new1)))

## Empty list should be accepted and leave things unchanged.
obj <- s$parameters
s$set_parameters(list())
expect_that(s$parameters, is_identical_to(obj))

## As should NULL, which is converted by R to list()
obj <- s$parameters
s$set_parameters(NULL)
expect_that(s$parameters, is_identical_to(obj))

ctrl <- new(Control)
expect_that(s$control$parameters, is_identical_to(ctrl$parameters))
ctrl$set_parameters(list(cohort_gradient_richardson=TRUE))
s$control <- ctrl
expect_that(s$control$parameters, is_identical_to(ctrl$parameters))

## Negative parameters should cause failure:
test_that("Negative parameters can't be set", {
  new.err <- list(lma=-1e-7)
  expect_that(new(Strategy, new.err), throws_error())
  old <- s$parameters
  expect_that(s$set_parameters(new.err), throws_error())
  expect_that(s$parameters, is_identical_to(old))
})

test_that("Copying a Strategy works", {
  s1 <- new(Strategy)
  s1$set_parameters(list(lma=1))
  s2 <- s1$copy()
  s2$set_parameters(list(lma=2))

  s1.p <- s1$parameters
  s2.p <- s2$parameters
  expect_that(s1.p[names(s1.p) != "lma"],
              is_identical_to(s2.p[names(s2.p) != "lma"]))
  expect_that(s1.p$lma, equals(1))
  expect_that(s2.p$lma, equals(2))
})
