source("helper-tree.R")

context("CohortSchedule")

n.types <- 2
sched <- new(CohortSchedule, n.types)

expect_that(sched$size, equals(0))
expect_that(sched$types, equals(n.types))
## Check empty event:
e <- sched$next_event
expect_that(e$time, is_identical_to(Inf))
expect_that(e$cohort, is_identical_to(NA_integer_))

set.seed(1)
t1 <- runif(10)
t2 <- runif(12)

## Times not sorted
expect_that(sched$set_times(t1, 1L), throws_error())
t1 <- sort(t1)
t2 <- sort(t2)
## Index out of bounds
expect_that(sched$set_times(t1, 0), throws_error())
expect_that(sched$set_times(t1, n.types + 1), throws_error())

## This will work fine though.
sched$set_times(t1, 1)
expect_that(sched$size, equals(length(t1)))

e <- sched$next_event
expect_that(e$time, is_identical_to(t1[[1]]))
cmp <- c()
for (i in seq_len(length(t1) + 1)) {
  e <- sched$next_event
  sched$pop()
  cmp <- c(cmp, e$time)
}
expect_that(cmp, is_identical_to(c(t1, Inf)))

expect_that(sched$times(1), equals(t1))

sched$set_times(t2, 2)
expect_that(sched$times(1), equals(t1))
expect_that(sched$times(2), equals(t2))
expect_that(sched$size, equals(length(t1) + length(t2)))

## Come up with the big list (2nd argument ensures stable sort, I
## hope).
tmp <- rbind(data.frame(time=t1, cohort=1),
             data.frame(time=t2, cohort=2))
tmp <- tmp[order(tmp$t, -tmp$cohort),]
cmp.cohort <- integer(0)
cmp.time   <- numeric(0)
sched$reset()
for (i in seq_len(length(t1) + length(t2) + 1)) {
  e <- sched$next_event
  sched$pop()
  cmp.time <- c(cmp.time, e$time)
  cmp.cohort <- c(cmp.cohort, e$cohort)
}

expect_that(cmp.time, equals(c(tmp$time, Inf)))
expect_that(cmp.cohort, equals(c(tmp$cohort-1, NA)))

expect_that(sched$next_time, equals(Inf))
sched$reset()
expect_that(sched$next_time, equals(tmp$time[[1]]))

## Check that resettting the times replaces them...
sched$set_times(t1 * 2, 1)
expect_that(sched$times(1), equals(2*t1))

sched0 <- new(CohortSchedule, 2)
expect_that(sched0$next_time, equals(Inf))

sched0 <- new(CohortSchedule, 0)
expect_that(sched0$next_time, equals(Inf))
