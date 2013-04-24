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
new.p <- list(mean_disturbance_interval=20.0)
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

