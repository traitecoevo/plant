source("helper-tree.R")

context("EBT")

p <- new(Parameters)
p$add_strategy(new(Strategy))

ebt <- new(EBT, p)

expect_that(inherits(ebt$patch, "Rcpp_PatchCohortTop"), is_true())

r <- pi/2
ebt$seed_rain <- seed_rain(r)

sched <- ebt$cohort_schedule
test_that("Empty CohortSchedule", {
  expect_that(sched$size, equals(0))
  expect_that(sched$next_time, equals(Inf))
})

## If the schedule is for the wrong number of species, it should cause
## an error...
expect_that(ebt$cohort_schedule <- new(CohortSchedule, sched$n.species + 1),
            throws_error())

t <- seq(0, 5, length=11)
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

ebt$reset()
test_that("EBT reset successful", {
  expect_that(ebt$time, equals(0))
  expect_that(ebt$patch$time, equals(0))
  expect_that(ebt$patch$ode_size, equals(0))
  expect_that(ebt$patch$n_individuals, equals(0))
})
ebt$seed_rain <- seed_rain(r)

## Run the whole schedule, and compare with a Patch run manually.
## TODO: This is really very slow at the moment...
patch <- new(PatchCohortTop, p)
patch$seed_rain <- ebt$seed_rain
species.index <- 1
sched$reset()
## Advance to time 't', then add species.
tt.p <- hh.p <- NULL
while (sched$remaining > 0) {
  patch$run_deterministic(sched$next_time)
  patch$add_seedling(species.index)
  sched$pop()
  tt.p <- c(tt.p, patch$time)
  hh.p <- c(hh.p, list(patch$height[[species.index]]))
}
n <- length(hh.p[[length(hh.p)]])
hh.p <- t(sapply(hh.p, function(x) c(x, rep(NA, n-length(x)))))

## Currently reset
## ebt <- new(EBT, p)
## ebt$seed_rain <- seed_rain(r)
## sched$reset()
## ebt$cohort_schedule <- sched
ebt$reset()
ebt$seed_rain <- seed_rain(r)
tt.e <- hh.e <- NULL
repeat {
  if (inherits(try(ebt$run_next(), silent=TRUE), "try-error"))
    break
  tt.e <- c(tt.e, ebt$time)
  hh.e <- c(hh.e, list(ebt$patch$height[[species.index]]))
}
n <- length(hh.e[[length(hh.e)]])
hh.e <- t(sapply(hh.e, function(x) c(x, rep(NA, n-length(x)))))

## However, *this* works:
ebt$reset()
ebt$seed_rain <- ebt$seed_rain
tt.e2 <- hh.e2 <- NULL
repeat {
  if (inherits(try(ebt$run_next(), silent=TRUE), "try-error"))
    break
  tt.e2 <- c(tt.e2, ebt$time)
  hh.e2 <- c(hh.e2, list(ebt$patch$height[[species.index]]))
}
n <- length(hh.e2[[length(hh.e2)]])
hh.e2 <- t(sapply(hh.e2, function(x) c(x, rep(NA, n-length(x)))))

expect_that(hh.e, equals(hh.e2))
expect_that(hh.p, equals(hh.e2))

## So - unlikely to be caused by the seeds (which are regenerated) or
## the environment (which is regenerated).  Not sure.  *Could* be
## something to do with the step size calculations.  Or the estimated
## slope?
if (interactive()) {
  matplot(tt.p, hh.p, type="o", col="black", pch=1, lty=1)
  matpoints(tt.e, hh.e, col="red", cex=.5, pch=1, type="o", lty=1)
}

rm(ebt)
gc()
