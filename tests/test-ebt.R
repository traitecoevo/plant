source("helper-tree.R")

context("EBT")

p <- new(Parameters)
p$add_strategy(new(Strategy))

ebt <- new(EBT, p)

patch <- ebt$patch
expect_that(inherits(patch, "Rcpp_PatchCohortTop"), is_true())

sched <- ebt$cohort_schedule
expect_that(sched$size, equals(0))
expect_that(sched$next_time, equals(Inf))

rm(ebt)
gc()
