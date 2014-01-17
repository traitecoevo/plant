source("helper-tree.R")

context("CohortSchedule")

n.species <- 2
sched <- new(CohortSchedule, n.species)

test_that("Empty CohortSchedule looks correct", {
  expect_that(sched$size,        equals(0))
  expect_that(sched$n_species,   equals(n.species))
  expect_that(sched$remaining,   equals(0))
  expect_that(sched$max_time,    equals(Inf))
  expect_that(sched$next_event,  throws_error())
  expect_that(sched$fixed_times, is_false())
  expect_that(sched$ode_times,   equals(numeric(0)))
})

set.seed(1)
t1 <- c(0.0, runif(10))
t2 <- c(0.0, runif(12))

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
  species.index <- 1
  expect_that(sched$times(species.index), equals(t1))
  expect_that(sched$fixed_times, is_false())
  e <- sched$next_event
  expect_that(e$time_introduction,   is_identical_to(t1[[1]]))
  expect_that(e$species_index,       equals(species.index))
})

drain.schedule <- function(sched) {
  sched$reset()
  cmp <- vector("list", sched$size)
  for (i in seq_len(sched$size)) {
    e <- sched$next_event
    cmp[[i]] <- c(e$species_index,
                  e$time_introduction,
                  e$times,
                  e$time_end)
    sched$pop()
  }
  if (!sched$fixed_times) {
    if (!all(sapply(cmp, length) == 5))
      stop("Expected exactly five elements for each schedule")
    cmp <- do.call(rbind, cmp)
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
max.t <- max.t <- max(c(t1, t2)) + mean(diff(sort(c(t1, t2))))
sched$max_time <- max.t

## Come up with the expected times (2nd argument ensures stable sort)
expected <- rbind(data.frame(species_index=1, start=t1),
                  data.frame(species_index=2, start=t2))
expected <- expected[order(expected$start, -expected$species_index),]
expected$end <- c(expected$start[-1], max.t)

test_that("Schedule contains correct information (two species)", {
  cmp <- drain.schedule(sched)

  expect_that(cmp[,1], equals(expected$species_index))
  expect_that(cmp[,2], equals(expected$start))
  expect_that(cmp[,3], equals(expected$start))
  expect_that(cmp[,4], equals(expected$end))
  expect_that(cmp[,5], equals(expected$end))
})

## Check that we behave sensibly at the end:
test_that("Schedule ends gracefully", {
  expect_that(sched$next_event, throws_error())
  expect_that(sched$max_time, equals(max.t))
  sched$reset()
  expect_that(sched$next_event$time_introduction,
              is_identical_to(min(c(t1, t2))))
})

test_that("Resetting the times replaces them", {
  t1.new <- t1 * .9
  sched$set_times(t1.new, 1)
  expect_that(sched$times(1), equals(t1.new))
  sched$set_times(t1, 1)
  expect_that(sched$times(1), equals(t1))
})

test_that("Setting max time behaves sensibly", {
  sched <- new(CohortSchedule, 2)
  sched$set_times(t1, 1)

  last.event <- function(x) {
    x <- x$copy()
    while (x$remaining > 1)
      x$pop()
    x$next_event
  }

  ## Before setting max_time, the finishing time will be Inf:
  e <- last.event(sched)
  expect_that(e$time_introduction, equals(last(t1)))
  expect_that(e$time_end,          equals(Inf))

  ## Set max_time to something stupid:
  expect_that(sched$max_time <- 0.5,            throws_error())
  expect_that(sched$max_time <- max(t1) - 1e-8, throws_error())

  ## And to something sensible:
  max.t <- max(t1) + 0.1
  sched$max_time <- max.t

  ## Make sure that the last event has been modified:
  e <- last.event(sched)
  expect_that(e$time_introduction, equals(last(t1)))
  expect_that(e$time_end,          equals(max.t))

  ## Now this will fail
  expect_that(sched$set_times(t1 * 2, 1), throws_error())
})

test_that("Bulk get/set of times works", {
  n <- 3
  sched <- new(CohortSchedule, n)

  set.seed(1)
  t.new <- lapply(seq_len(n), function(...) sort(runif(rpois(1, 10))))
  sched$all_times <- t.new
  expect_that(sched$all_times, is_identical_to(t.new))

  expect_that(sched$all_times <- t.new[1], throws_error())
  expect_that(sched$all_times <- t.new[[1]], throws_error())
  expect_that(sched$all_times <- t.new[c(1, seq_len(n))], throws_error())
  ## Unfortunately, this does not throw.
  ## expect_that(sched$all_times <- seq_len(n), throws_error())
})

## Fixed times.  First set some impossible cases:
test_that("Malformed ode_times objects are rejected", {
  ## Too few values:
  expect_that(sched$ode_times <- numeric(0), throws_error())
  expect_that(sched$ode_times <- c(0.0),     throws_error())
  ## Does not start at 0
  expect_that(sched$ode_times <- c(1, 2, 3), throws_error())
  ## Does not finish at time_max
  expect_that(sched$ode_times <- c(0.0, 2, 3), throws_error())
  ## Is not sorted:
  expect_that(sched$ode_times <- sched$max_time * c(0, .5, .3, 1),
              throws_error())
  ## ...and check that none of these caused the times to be set
  expect_that(sched$fixed_times, is_false())
  expect_that(sched$ode_times,   equals(numeric(0)))
})

## So, now, manually get the times set up.  For a real challenge, we
## should add some more exact hits in here.  Building the expected
## times in R is hard enough!
t.ode <- seq(0, max.t, length=14)
idx <- findInterval(t.ode, c(expected$start, max.t), TRUE)
tmp <- unname(t(apply(expected[c("start", "end")], 1, unlist)))

expected.ode <- lapply(seq_len(nrow(expected)), function(i)
                       c(tmp[i,1],
                         setdiff(t.ode[idx == i], tmp[i,]),
                         tmp[i,2]))

## New schedule because setting and resetting may have changed cohort
## order.
sched <- new(CohortSchedule, n.species)
sched$set_times(t1, 1L)
sched$set_times(t2, 2L)
sched$max_time <- max.t
sched$ode_times <- t.ode

test_that("Times were set correctly", {
  expect_that(sched$fixed_times, is_true())
  expect_that(sched$ode_times,   is_identical_to(t.ode))

  cmp <- drain.schedule(sched)

  expect_that(sapply(cmp, first),  equals(expected$species_index))
  expect_that(sapply(cmp, second), equals(expected$start))
  expect_that(sapply(cmp, last),   equals(expected$end))

  expect_that(lapply(cmp, function(x) x[3:(length(x) - 1)]),
              equals(expected.ode))
})

## check we can clear times:
sched$clear_ode_times()
test_that("Times were cleared", {
  expect_that(sched$fixed_times, is_false())
  expect_that(sched$ode_times,   equals(numeric(0)))
})

sched$max_time <- Inf
sched$ode_times <- t.ode
test_that("Times were set correctly", {
  expect_that(sched$fixed_times, is_true())
  expect_that(sched$ode_times,   is_identical_to(t.ode))
  expect_that(sched$max_time,    is_identical_to(max(t.ode)))
})

test_that("State get/set works correctly", {
  sched$reset()
  state <- sched$state
  # lapply(seq_len(sched$n_species), function(i) sched$times(i)))
  expect_that(state$times, is_identical_to(list(t1, t2)))
  expect_that(state$ode_times, is_identical_to(t.ode))
  expect_that(state$remaining, equals(sched$size))

  for (i in 1:10)
    sched$pop()
  state <- sched$state
  expect_that(state$times, is_identical_to(list(t1, t2)))
  expect_that(state$ode_times, is_identical_to(t.ode))
  expect_that(state$remaining, equals(sched$remaining))

  sched2 <- new(CohortSchedule, n.species)
  sched2$state <- state
  expect_that(sched2$state, is_identical_to(state))
  expect_that(sched2$remaining, is_identical_to(sched$remaining))

  ## Repeat without ode times:
  sched$clear_ode_times()
  sched$reset()
  state <- sched$state
  expect_that(state$times, is_identical_to(list(t1, t2)))
  expect_that(state$ode_times, equals(numeric(0)))
  expect_that(state$remaining, equals(sched$size))

  for (i in 1:10)
    sched$pop()
  state <- sched$state
  expect_that(state$times, is_identical_to(list(t1, t2)))
  expect_that(state$ode_times, is_identical_to(numeric(0)))
  expect_that(state$remaining, equals(sched$remaining))

  sched2 <- new(CohortSchedule, n.species)
  sched2$state <- state
  expect_that(sched2$state, is_identical_to(state))
  expect_that(sched2$remaining, is_identical_to(sched$remaining))
})

test_that("Can expand CohortSchedule", {
  sched <- new(CohortSchedule, 1)
  max.t <- 10
  times1 <- sort(runif(10))
  sched$max_time <- max.t
  sched$set_times(times1, 1)

  expect_that(sched$n_species, equals(1))
  expect_that(sched$max_time, is_identical_to(max.t))
  expect_that(sched$times(1), is_identical_to(times1))

  times2 <- sort(runif(20))
  sched3 <- sched$expand(2, times2) # expand by two species:
  expect_that(sched3$n_species, equals(3))
  expect_that(sched3$max_time, is_identical_to(max.t))
  expect_that(sched3$times(1), is_identical_to(times1))
  expect_that(sched3$times(2), is_identical_to(times2))
  expect_that(sched3$times(3), is_identical_to(times2))

  # Schedules are independent:
  sched3$max_time <- 2 * max.t
  expect_that(sched$max_time,  is_identical_to(max.t))
  expect_that(sched3$max_time, is_identical_to(2 * max.t))
})
