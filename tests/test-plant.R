source("helper-tree.R")

context("Plant")

s <- new(Strategy)

## Create a plant, and delete it.  Should not crash.
p <- new(Plant, s)
rm(p)
gc()

