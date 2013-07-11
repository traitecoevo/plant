source("helper-tree.R")

context("CohortSchedule")

n.species <- 2
sched <- new(CohortSchedule, n.species)

test_that("Empty CohortSchedule looks correct", {
  expect_that(sched$size, equals(0))
  expect_that(sched$n_species, equals(n.species))
  expect_that(sched$remaining, equals(0))
  expect_that(sched$max_time, equals(Inf))
})

## Check empty event:
test_that("Empty event looks correct", {
  e <- sched$next_event
  expect_that(e$time, is_identical_to(Inf))
  expect_that(e$cohort, is_identical_to(NA_integer_))
})

set.seed(1)
t1 <- runif(10)
t2 <- runif(12)

## Times not sorted
test_that("Unsorted times cause error", {
  expect_that(sched$set_times(t1, 1L), throws_error())
})

t1 <- sort(t1)
t2 <- sort(t2)

test_that("Negative times cause error", {
  expect_that(sched$set_times(c(-.1, t1), 1L), throws_error())
})

test_that("Index out of bounds throws error", {
  expect_that(sched$set_times(t1, 0), throws_error())
  expect_that(sched$set_times(t1, n.species + 1), throws_error())
})

sched$set_times(t1, 1)
test_that("Can set cohort times", {
  expect_that(sched$size, equals(length(t1)))
  expect_that(sched$remaining, equals(length(t1)))
})

e <- sched$next_event
expect_that(e$time, is_identical_to(t1[[1]]))

cmp <- c()
sched$reset()
for (i in seq_len(length(t1) + 1)) {
  e <- sched$next_event
  cmp <- c(cmp, e$time)
  if (sched$remaining > 0)
    sched$pop()
}
expect_that(cmp, is_identical_to(c(t1, Inf)))
expect_that(sched$remaining, equals(0))

expect_that(sched$times(1), equals(t1))

sched$set_times(t2, 2)
expect_that(sched$times(1), equals(t1))
expect_that(sched$times(2), equals(t2))
expect_that(sched$size, equals(length(t1) + length(t2)))

## Force a max_time for the next run through:
max.t <- max(c(t1, t2) + 2)
sched$max_time <- max.t

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
  cmp.time <- c(cmp.time, e$time)
  cmp.cohort <- c(cmp.cohort, e$cohort)
  if (sched$remaining > 0)
    sched$pop()
}

expect_that(cmp.time, equals(c(tmp$time, max.t)))
expect_that(cmp.cohort, equals(c(tmp$cohort, NA)))

expect_that(sched$next_time, equals(max.t))
sched$reset()
expect_that(sched$next_time, equals(tmp$time[[1]]))

drop <- 3
for (i in seq_len(drop))
  sched$pop()
expect_that(sched$remaining, equals(length(t1) + length(t2) - drop))

## Check that resettting the times replaces them...
sched$set_times(t1 * 2, 1)
expect_that(sched$times(1), equals(2*t1))

sched0 <- new(CohortSchedule, 2)
expect_that(sched0$next_time, equals(Inf))

sched0 <- new(CohortSchedule, 0)
expect_that(sched0$next_time, equals(Inf))

## Give the logic around max_t a better workout:
sched <- new(CohortSchedule, 2)
sched$set_times(t1, 1)
expect_that(sched$max_time <- 0.5, throws_error())
expect_that(sched$max_time <- max(t1) - 1e-8, throws_error())
sched$max_time <- max(t1)
## Now this will fail
expect_that(sched$set_times(t1 * 2, 1), throws_error())
