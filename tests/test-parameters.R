source("helper-tree.R")

context("Parameters")

p <- new(Parameters)

## Must start empty
expect_that(p$size, is_identical_to(0L))

## And these are the defaults:
expected <- list(c_ext=0.5,
                 mean_disturbance_interval=30.0)
expect_that(p$get_parameters(), is_identical_to(expected))

## Set a parameter and check that it is actually set
new.p <- list(mean_disturbance_interval=20.0)
expected.s <- modifyList(expected, new.p)
p$set_parameters(new.p)
expect_that(p$get_parameters(), is_identical_to(expected.s))

## Set an invalid key and check that it throws error:
new.err <- list(unknown_key=1)
expect_that(p$set_parameters(new.err), throws_error())
expect_that(p$get_parameters(), is_identical_to(expected.s))

## Getting a nonexistant strategy should cause an error
expect_that(p$get_strategy(0L), throws_error())
expect_that(p$get_strategy(10), throws_error())

## Add a strategy:
p$add_strategy(list())

expect_that(p$size, is_identical_to(1L))
## Should not have changed any parameters
expect_that(p$get_parameters(), is_identical_to(expected.s))

## The added strategy should be the same as the default strategy:
s <- new(Strategy)
cmp <- s$get_parameters()
expect_that(p$get_strategy(0L), is_identical_to(cmp))
expect_that(p$get_strategies(), is_identical_to(list(cmp)))

new.s <- list(c_acc=4.1, lma=1.2)
p$set_strategy(new.s, 0L)

s$set_parameters(new.s)
cmp.mod <- s$get_parameters()
expect_that(p$get_strategy(0L), is_identical_to(cmp.mod))
expect_that(p$get_strategies(), is_identical_to(list(cmp.mod)))

## Add another strategy
p$add_strategy(list())
expect_that(p$size, is_identical_to(2L))
expect_that(p$get_strategy(1L), is_identical_to(cmp))
expect_that(p$get_strategies(), is_identical_to(list(cmp.mod, cmp)))

expect_that(p$set_strategy(new.s, 10L),
            throws_error())
