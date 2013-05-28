source("helper-tree.R")

context("EBT")

p <- new(Parameters)
p$add_strategy(new(Strategy))

ebt <- new(EBT, p)

patch <- ebt$patch
expect_that(inherits(patch, "Rcpp_PatchCohortTop"), is_true())

rm(ebt)
gc()
