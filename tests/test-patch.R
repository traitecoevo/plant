source("helper-tree.R")

context("Patch")

p <- new(Parameters)
p$add_strategy(list())

## TODO: This should probably give back the C++ object?
s.p <- p$get_strategy(0)
s <- new(Strategy)
s$set_parameters(s.p)

popn <- new(Patch, p)

expect_that(popn$size, equals(p$size))

## We've not added any seeds yet, so this should be the empty list
plants <- popn$get_plants(0L)
expect_that(plants, is_identical_to(list()))

popn$add_seed(0)

plants <- popn$get_plants(0L)
expect_that(length(plants), equals(1))
expect_that(length(plants[[1]]), equals(1))

cmp <- new(Plant, s)

expect_that(plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

rm(popn)
gc()
