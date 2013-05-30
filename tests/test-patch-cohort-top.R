source("helper-tree.R")

context("Patch [CohortTop]")

p <- new(Parameters)
p$add_strategy(new(Strategy))

## A plant that will be the same in terms of strategy (and initial
## mass).
cmp <- new(Plant, p[[1]])

patch.p <- new(Patch,  p)
patch.c <- new(PatchCohortTop, p)

patch.c$add_seedling(1)
patch.c$ode_size

y <- patch.c$ode_values
## This gives a "requested rain out of bounds", which is probably the
## cause of the EBT error -- nice.
##
## So, update to get the seed rain initialised, and we're probably
## sorted.
##
## I'm still a bit confused though, as to why this crashes the ebt...
## patch.c$set_ode_values(0, y)
