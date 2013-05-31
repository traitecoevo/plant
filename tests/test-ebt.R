source("helper-tree.R")

context("EBT")

p <- new(Parameters)
p$add_strategy(new(Strategy))

ebt <- new(EBT, p)

expect_that(inherits(ebt$patch, "Rcpp_PatchCohortTop"), is_true())

r <- pi/2
ebt$set_seed_rain(seed_rain(r))

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

## This will fail:
ebt$patch$n_individuals
ebt$run_next()
expect_that(ebt$cohort_schedule$remaining, equals(length(t) - 1))
## This should be 4
ebt$patch$ode_size

## This will now fail:
## ebt$run_next()
## Failure is in CohortTop::set_ode_values still.

rm(ebt)
gc()
