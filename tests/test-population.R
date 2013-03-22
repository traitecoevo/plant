source("helper-tree.R")

context("Population")

p <- new(Parameters)
p$add_strategy(list())

popn <- new(Population, p)

rm(popn)
gc()
