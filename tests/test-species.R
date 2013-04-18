source("helper-tree.R")
options(error=traceback)

context("Species [Plant]")

s <- new(Strategy)

sp <- new(Species, s)
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))

sp$add_seeds(1)
expect_that(sp$size, equals(1))
expect_that(sp$n_individuals, equals(1))

rm(sp)
gc()

