source("helper-tree.R")

context("Patch")

p <- new(Parameters)
p$add_strategy(list())

popn <- new(Patch, p)

rm(popn)
gc()
