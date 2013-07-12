source("helper-tree.R")

context("CohortSchedule")

n.species <- 2
sched <- new(CohortSchedule, n.species)

test_that("Empty CohortSchedule looks correct", {
  expect_that(sched$size, equals(0))
  expect_that(sched$n_species, equals(n.species))
  expect_that(sched$remaining, equals(0))
  expect_that(sched$max_time, equals(Inf))
  expect_that(sched$next_event, throws_error())
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

test_that("Schedule looks correctly set up", {
  expect_that(sched$times(1), equals(t1))
  expect_that(sched$next_event$time_introduction,
              is_identical_to(t1[[1]]))
})

drain.schedule <- function(sched) {
  sched$reset()
  cmp <- matrix(nrow=sched$size, ncol=5)
  for (i in seq_len(sched$size)) {
    e <- sched$next_event
    cmp[i,] <- c(e$species_index,
                 e$time_introduction,
                 e$times,
                 e$time_end)
    sched$pop()
  }
  cmp
}

test_that("Schedule contains correct information (one species)", {
  cmp <- drain.schedule(sched)
  expect_that(sched$remaining, equals(0))
  n <- length(t1)
  expect_that(cmp[,1], equals(rep(1, n)))
  expect_that(cmp[,2], equals(t1))
  expect_that(cmp[,3], equals(t1))
  expect_that(cmp[,4], equals(c(t1[-1], Inf)))
  expect_that(cmp[,5], equals(c(t1[-1], Inf)))
})

## Now, set times for the second species:
sched$set_times(t2, 2)
test_that("Schedule correctly set up with two species", {
  expect_that(sched$times(1), equals(t1))
  expect_that(sched$times(2), equals(t2))
  expect_that(sched$size, equals(length(t1) + length(t2)))
})

## Force a max_time for this run through:
max.t <- max(c(t1, t2) + 2)
sched$max_time <- max.t

test_that("Schedule contains correct information (two species)", {
  cmp <- drain.schedule(sched)

  ## Come up with the big list (2nd argument ensures stable sort)
  tmp <- rbind(data.frame(species_index=1, time=t1),
               data.frame(species_index=2, time=t2))
  tmp <- tmp[order(tmp$t, -tmp$species_index),]

  expect_that(cmp[,1], equals(tmp$species_index))
  expect_that(cmp[,2], equals(tmp$time))
  expect_that(cmp[,3], equals(tmp$time))
  expect_that(cmp[,4], equals(c(tmp$time[-1], max.t)))
  expect_that(cmp[,5], equals(c(tmp$time[-1], max.t)))
})

## Check that we behave sensibly at the end:
test_that("Schedule ends gracefully", {
  expect_that(sched$next_event, throws_error())
  expect_that(sched$max_time, equals(max.t))
  sched$reset()
  expect_that(sched$next_event$time_introduction,
              is_identical_to(min(c(t1, t2))))
})

test_that("Resettting the times replaces them", {
  sched$set_times(t1 * 2, 1)
  expect_that(sched$times(1), equals(2*t1))
})

test_that("Setting max time behaves sensibly", {
  sched <- new(CohortSchedule, 2)
  sched$set_times(t1, 1)
  expect_that(sched$max_time <- 0.5, throws_error())
  expect_that(sched$max_time <- max(t1) - 1e-8, throws_error())
  sched$max_time <- max(t1)
  ## Now this will fail
  expect_that(sched$set_times(t1 * 2, 1), throws_error())
})
