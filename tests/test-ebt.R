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

## If the schedule is for the wrong number of species, it should cause
## an error...
expect_that(ebt$cohort_schedule <- new(CohortSchedule, sched$n.species + 1),
            throws_error())

t <- seq(0, 10, length=11)
sched$set_times(t, 1)
ebt$cohort_schedule <- sched

expect_that(ebt$time, equals(0.0))

rm(ebt)
gc()
