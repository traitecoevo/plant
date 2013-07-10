source("helper-tree.R")

context("EBT")

p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- pi/2

ebt <- new(EBT, p)

## Check that the underlying Patch really is a Patch<CohortTop>, not
## one of the other Patch types (which would just fail miserably below
## here, if they even compile).
expect_that(inherits(ebt$patch, "Rcpp_PatchCohortTop"), is_true())

sched <- ebt$cohort_schedule
test_that("Empty CohortSchedule", {
  expect_that(sched$size, equals(0))
  expect_that(sched$next_time, equals(Inf))
})

## If the schedule is for the wrong number of species, it should cause
## an error...
sched2 <- new(CohortSchedule, sched$n_species + 1)
expect_that(ebt$cohort_schedule <- sched2, throws_error())

## Build a schedule for 11 introductions from t=0 to t=5
t <- seq(0, 5, length=11)
sched$set_times(t, 1)
ebt$cohort_schedule <- sched

## Before starting, check that the EBT is actually empty
test_that("EBT starts empty", {
  expect_that(ebt$time, equals(0.0))
  expect_that(ebt$patch$ode_size, equals(0))
})

## Check that we can advance through and add the cohort at time zero.
ebt$patch$n_individuals
ebt$run_next()

test_that("EBT adds cohort successfully", {
  expect_that(ebt$cohort_schedule$remaining, equals(length(t) - 1))
  expect_that(ebt$time, equals(t[1]))
  expect_that(ebt$patch$ode_size, equals(4))
})

## Continue up to the introduction of cohort 2 (this will actually
## advance time)
ebt$run_next()
test_that("EBT ran successfully", {
  expect_that(ebt$cohort_schedule$remaining, equals(length(t) - 2))
  expect_that(ebt$time, equals(t[2]))
  expect_that(ebt$patch$ode_size, equals(4*2))
})

## Reset everything
ebt$reset()
test_that("EBT reset successful", {
  expect_that(ebt$time, equals(0))
  expect_that(ebt$patch$time, equals(0))
  expect_that(ebt$patch$ode_size, equals(0))
  expect_that(ebt$patch$n_individuals, equals(0))
})

## Run the whole schedule using a Patch<CohortTop>, manually moving
## things along the schedule.
patch <- new(PatchCohortTop, p)
species.index <- 1
sched$reset()
## Advance to time 't', then add species.
tt.p <- hh.p <- NULL
solver <- solver.from.odetarget(patch, p$control$ode_control)
while (sched$remaining > 0) {
  solver$advance(sched$next_time)
  patch$add_seedling(species.index)
  solver$set_state(patch$ode_values, patch$time)
  sched$pop()
  tt.p <- c(tt.p, patch$time)
  hh.p <- c(hh.p, list(patch$height[[species.index]]))
}
n <- length(hh.p[[length(hh.p)]])
hh.p <- t(sapply(hh.p, function(x) c(x, rep(NA, n-length(x)))))

test_that("Run looks successful", {
  expect_that(n, equals(length(t)))
  expect_that(all(diag(hh.p) == new(Plant, p[[1]])$height), is_true())
})

## Run the whole schedule using the EBT.
ebt$reset()
tt.e <- hh.e <- NULL
repeat {
  if (inherits(try(ebt$run_next(), silent=TRUE), "try-error"))
    break
  tt.e <- c(tt.e, ebt$time)
  hh.e <- c(hh.e, list(ebt$patch$height[[species.index]]))
}
n <- length(hh.e[[length(hh.e)]])
hh.e <- t(sapply(hh.e, function(x) c(x, rep(NA, n-length(x)))))

expect_that(hh.e, is_identical_to(hh.e))

## Then, check that resetting the cohort allows rerunning easily:
ebt$reset()
tt.e2 <- hh.e2 <- NULL
repeat {
  if (inherits(try(ebt$run_next(), silent=TRUE), "try-error"))
    break
  tt.e2 <- c(tt.e2, ebt$time)
  hh.e2 <- c(hh.e2, list(ebt$patch$height[[species.index]]))
}
n <- length(hh.e2[[length(hh.e2)]])
hh.e2 <- t(sapply(hh.e2, function(x) c(x, rep(NA, n-length(x)))))

expect_that(hh.e, is_identical_to(hh.e2))

if (interactive()) {
  matplot(tt.p, hh.p, type="o", col="black", pch=1, lty=1)
  matpoints(tt.e, hh.e, col="red", cex=.5, pch=1, type="o", lty=1)
}

rm(ebt)
gc()
