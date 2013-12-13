source("helper-tree.R")

context("Parameters")

p <- new(Parameters)

## Must start empty
expect_that(p$size, equals(0))

## And these are the defaults:
expected <- list(Pi_0=0.25,
                 c_ext=0.5,
                 mean_disturbance_interval=30.0,
                 n_patches=10,
                 patch_area=10.0)
expect_that(p$parameters, is_identical_to(expected))

## Set a parameter and check that it is actually set
new.p <- list(patch_area=20.0)
expected.s <- modifyList(expected, new.p)
p$set_parameters(new.p)
expect_that(p$parameters, is_identical_to(expected.s))

## Set an invalid key and check that it throws error:
new.err <- list(unknown_key=1)
expect_that(p$set_parameters(new.err), throws_error())
expect_that(p$parameters, is_identical_to(expected.s))

## Getting a nonexistant strategy should cause an error
expect_that(p[[1]],  throws_error())
expect_that(p[[10]], throws_error())

## Add a (default) strategy:
p$add_strategy(new(Strategy))

expect_that(p$size, equals(1))
## Should not have changed any parameters
expect_that(p$parameters, is_identical_to(expected.s))

## The added strategy should be the same as the default strategy:
s <- new(Strategy)
cmp <- s$parameters
expect_that(p[[1]]$parameters, is_identical_to(cmp))
res <- p$strategies
expect_that(length(res), equals(1))
expect_that(res[[1]]$parameters, is_identical_to(cmp))

mod <- list(c_acc=4.1, lma=1.2)
new.s <- new(Strategy, mod)
cmp.mod <- new.s$parameters
## Quick check
expect_that(cmp.mod,
            equals(modifyList(s$parameters, mod)))

p$add_strategy(new.s)

res <- p$strategies
expect_that(length(res), equals(2))
expect_that(res[[1]]$parameters, is_identical_to(cmp))
expect_that(res[[2]]$parameters, is_identical_to(cmp.mod))

expect_that(p[[1]]$parameters, is_identical_to(res[[1]]$parameters))
expect_that(p[[2]]$parameters, is_identical_to(res[[2]]$parameters))

ctrl <- new(Control)
expect_that(p$control$parameters, is_identical_to(ctrl$parameters))
ctrl$set_parameters(list(cohort_gradient_richardson=TRUE))
p$control <- ctrl
expect_that(p$control$parameters, is_identical_to(ctrl$parameters))

ctrl.extra <- list(plant_seed_tol=1)
p$set_control_parameters(ctrl.extra)
expect_that(p$control$parameters,
            is_identical_to(modifyList(ctrl$parameters, ctrl.extra)))
test_that("Control parameters propogate to strategies", {
  expect_that(p[[1]]$control$parameters[["plant_seed_tol"]],
              equals(ctrl.extra[["plant_seed_tol"]]))
})

test_that("Seed rain is correct length", {
  expect_that(length(p$seed_rain), equals(p$size))
  expect_that(p$seed_rain, is_identical_to(rep(1.0, p$size)))
})

test_that("Seed rain setting works as expected", {
  expect_that(p$seed_rain <- numeric(0), throws_error())
  expect_that(p$seed_rain <- 1,          throws_error())
  expect_that(p$seed_rain <- rep(1, 3),  throws_error())
  r <- runif(2)
  p$seed_rain <- r
  expect_that(p$seed_rain, is_identical_to(r))
})

test_that("Resident flag is correct length", {
  expect_that(length(p$is_resident), equals(p$size))
  expect_that(p$is_resident, is_identical_to(rep(TRUE, p$size)))
})

test_that("Resident flag setting works as expected", {
  expect_that(p$is_resident <- logical(0),   throws_error())
  expect_that(p$is_resident <- TRUE,         throws_error())
  expect_that(p$is_resident <- rep(TRUE, 3), throws_error())
  r <- c(FALSE, TRUE)
  p$is_resident <- r
  expect_that(p$is_resident, is_identical_to(r))
})

test_that("Can directly add mutant strategies", {
  r <- p$is_resident
  p$add_strategy_mutant(new(Strategy))
  expect_that(p$is_resident,
              is_identical_to(c(r, FALSE)))
})

test_that("Cloning a Parameters works", {
  p1 <- new(Parameters)
  p1$set_parameters(list(c_ext=1))
  p2 <- p1$clone()
  p2$set_parameters(list(c_ext=2))

  p1.p <- p1$parameters
  p2.p <- p2$parameters
  expect_that(p1.p[names(p1.p) != "c_ext"],
              is_identical_to(p2.p[names(p2.p) != "c_ext"]))
  expect_that(p1.p$c_ext, equals(1))
  expect_that(p2.p$c_ext, equals(2))
})
