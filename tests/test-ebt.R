source("helper-tree.R")

context("EBT")

p <- new(Parameters)
p$add_strategy(new(Strategy))

ebt <- new(EBT, p)

expect_that(inherits(ebt$patch, "Rcpp_PatchCohortTop"), is_true())

r <- pi/2
ebt$set_seed_rain(seed_rain(r))

sched <- ebt$cohort_schedule
test_that("Empty CohortSchedule", {
  expect_that(sched$size, equals(0))
  expect_that(sched$next_time, equals(Inf))
})

## If the schedule is for the wrong number of species, it should cause
## an error...
expect_that(ebt$cohort_schedule <- new(CohortSchedule, sched$n.species + 1),
            throws_error())

t <- seq(0, 10, length=11)
sched$set_times(t, 1)
ebt$cohort_schedule <- sched

test_that("EBT starts empty", {
  expect_that(ebt$time, equals(0.0))
  expect_that(ebt$patch$ode_size, equals(0))
})

ebt$patch$n_individuals
ebt$run_next()

test_that("EBT adds cohort successfully", {
  expect_that(ebt$cohort_schedule$remaining, equals(length(t) - 1))
  expect_that(ebt$time, equals(0.0))
  expect_that(ebt$patch$ode_size, equals(4))
})

## Run up to the introduction of cohort 2:
ebt$run_next()

test_that("EBT ran successfully", {
  expect_that(ebt$cohort_schedule$remaining, equals(length(t) - 2))
  expect_that(ebt$time, equals(t[2]))
  expect_that(ebt$patch$ode_size, equals(4*2))
})

rm(ebt)
gc()
