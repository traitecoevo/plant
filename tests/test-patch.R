source("helper-tree.R")

context("Patch")

p <- new(Parameters)
p$add_strategy(list())

## A plant that will be the same in terms of strategy (and initial
## mass).
cmp <- new(Plant, p$get_strategy(0))

popn <- new(Patch, p)

expect_that(popn$size, equals(p$size))

## We've not added any seeds yet, so this should be the empty list
plants <- popn$get_plants(0L)
expect_that(plants, is_identical_to(list()))

## And the height must be zero
expect_that(popn$height_max, is_identical_to(0.0))

## Add a single seed to this
popn$add_seed(0)

plants <- popn$get_plants(0L)
expect_that(length(plants), equals(1))
expect_that(length(plants[[1]]), equals(1))

expect_that(plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

expect_that(popn$height_max,
            is_identical_to(cmp$vars_size[["height"]]))

rm(popn)
gc()
