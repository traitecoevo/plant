drain_schedule <- function(sched) {
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
  if (!sched$use_ode_times) {
    if (!all(sapply(cmp, length) == 5)) {
      stop("Expected exactly five elements for each schedule")
    }
    cmp <- rbind_list(cmp)
  }
  cmp
}

context("CohortSchedule")

test_that("CohortScheduleEvent", {
  e <- CohortScheduleEvent(pi, 1)

  expect_identical(e$species_index, 1L)
  expect_identical(e$species_index_raw, 0.0)

  e$species_index <- 2L
  expect_identical(e$species_index, 2L)
  expect_identical(e$species_index_raw, 1.0)

  expect_identical(e$times, pi)
  expect_identical(e$time_introduction, pi)
  expect_identical(e$time_end, pi)
})

test_that("Empty CohortSchedule", {
  n_species <- 2
  sched <- CohortSchedule(n_species)

  expect_is(sched, "CohortSchedule")
  expect_equal(sched$size, 0)
  expect_equal(sched$n_species, n_species)

  expect_equal(sched$remaining, 0)
  expect_equal(sched$max_time, Inf)
  expect_error(sched$next_event)
  expect_false(sched$use_ode_times)
  expect_equal(sched$ode_times, numeric(0))
})

test_that("Corner cases", {
  n_species <- 2
  sched <- CohortSchedule(n_species)

  set.seed(1)
  t1 <- c(0.0, runif(10))
  t2 <- c(0.0, runif(12))

  expect_error(sched$set_times(t1, 1L), "Times must be sorted")

  t1 <- sort(t1)
  t2 <- sort(t2)

  expect_error(sched$set_times(t1, -1), "Invalid value for index")
  expect_error(sched$set_times(t1, 0), "Invalid value for index")
  expect_error(sched$set_times(t1, n_species + 1L), "Index 3 out of bounds")
})

test_that("Set times (one species)", {
  set.seed(1)
  t1 <- sort(c(0.0, runif(10)))
  t2 <- sort(c(0.0, runif(12)))

  n_species <- 2
  sched <- CohortSchedule(n_species)
  sched$set_times(t1, 1L)

  expect_equal(sched$size, length(t1))
  expect_equal(sched$remaining, length(t1))

  species_index <- 1
  expect_equal(sched$times(species_index), t1)
  expect_false(sched$use_ode_times)
  e <- sched$next_event
  expect_identical(e$time_introduction, t1[[1]])
  expect_equal(e$species_index, species_index)

  cmp <- drain_schedule(sched)

  expect_equal(sched$remaining, 0)
  n <- length(t1)
  expect_equal(cmp[,1], rep(1, n))
  expect_equal(cmp[,2], t1)
  expect_equal(cmp[,3], t1)
  expect_equal(cmp[,4], c(t1[-1], Inf))
  expect_equal(cmp[,5], c(t1[-1], Inf))
})

test_that("Set times (two species)", {
  set.seed(1)
  t1 <- sort(c(0.0, runif(10)))
  t2 <- sort(c(0.0, runif(12)))

  n_species <- 2
  sched <- CohortSchedule(n_species)
  sched$set_times(t1, 1L)
  sched$set_times(t2, 2L)

  expect_equal(sched$times(1), t1)
  expect_equal(sched$times(2), t2)
  expect_equal(sched$size, length(t1) + length(t2))

  ## Force a max_time for this run through:
  max_t <- max(c(t1, t2)) + mean(diff(sort(c(t1, t2))))
  sched$max_time <- max_t

  ## Come up with the expected times (2nd argument ensures stable sort)
  expected <- rbind(data.frame(species_index=1, start=t1),
                    data.frame(species_index=2, start=t2))
  expected <- expected[order(expected$start, -expected$species_index),]
  expected$end <- c(expected$start[-1], max_t)

  cmp <- drain_schedule(sched)

  expect_equal(cmp[,1], expected$species_index)
  expect_equal(cmp[,2], expected$start)
  expect_equal(cmp[,3], expected$start)
  expect_equal(cmp[,4], expected$end)
  expect_equal(cmp[,5], expected$end)

  expect_error(sched$next_event, "All events completed")
  expect_equal(sched$max_time, max_t)
  sched$reset()
  expect_identical(sched$next_event$time_introduction, min(c(t1, t2)))
})

test_that("Resetting times replaces them", {
  set.seed(1)
  t1 <- sort(c(0.0, runif(10)))
  t2 <- sort(c(0.0, runif(12)))

  n_species <- 2
  sched <- CohortSchedule(n_species)
  sched$set_times(t1, 1L)
  sched$set_times(t2, 2L)

  t1_new <- t1 * .9
  sched$set_times(t1_new, 1)
  expect_equal(sched$times(1), t1_new)
  sched$set_times(t1, 1)
  expect_equal(sched$times(1), t1)
})

test_that("Setting max time behaves sensibly", {
  set.seed(1)
  t1 <- sort(c(0.0, runif(10)))
  t2 <- sort(c(0.0, runif(12)))

  sched <- CohortSchedule(2)
  sched$set_times(t1, 1)

  last_event <- function(x) {
    x <- x$copy()
    while (x$remaining > 1L) {
      x$pop()
    }
    x$next_event
  }

  ## Before setting max_time, the finishing time will be Inf:
  e <- last_event(sched)
  expect_equal(e$time_introduction, last(t1))
  expect_equal(e$time_end, Inf)

  ## Set max_time to something stupid:
  expect_error(sched$max_time <- 0.5, "max_time must be at least the final")
  expect_error(sched$max_time <- max(t1) - 1e-8, "max_time must be at least the final")

  ## And to something sensible:
  max_t <- max(t1) + 0.1
  sched$max_time <- max_t

  ## Make sure that the last event has been modified:
  e <- last_event(sched)
  expect_equal(e$time_introduction, last(t1))
  expect_equal(e$time_end, max_t)

  ## Now this will fail
  expect_error(sched$set_times(t1 * 2, 1), "Times cannot be greater than max_time")
})

test_that("Bulk get/set of times works", {
  n <- 3
  sched <- CohortSchedule(n)

  set.seed(1)
  t_new <- lapply(seq_len(n), function(...) sort(runif(rpois(1, 10))))
  sched$all_times <- t_new
  expect_identical(sched$all_times, t_new)

  expect_error(sched$all_times <- t_new[1], "Incorrect length input")
  expect_error(sched$all_times <- t_new[[1]], "Incorrect length input")
  expect_error(sched$all_times <- t_new[c(1, seq_len(n))], "Incorrect length input")
  ## Unfortunately, this does not throw.
  ## expect_error(sched$all_times <- seq_len(n))
})

## Fixed times.  First set some impossible cases:
test_that("ode_times", {
  set.seed(1)
  t1 <- sort(c(0.0, runif(10)))
  t2 <- sort(c(0.0, runif(12)))
  n_species <- 2
  max_t <- max(c(t1, t2)) + mean(diff(sort(c(t1, t2))))

  sched <- CohortSchedule(n_species)
  sched$set_times(t1, 1L)
  sched$set_times(t2, 2L)
  sched$max_time <- max_t

  sched$ode_times <- numeric(0)
  expect_identical(sched$ode_times, numeric(0))

  ## Too few values:
  expect_error(sched$ode_times <- c(0.0), "Need at least two times")
  ## Does not start at 0
  expect_error(sched$ode_times <- c(1, 2, 3), "First time must be exactly zero")
  ## Does not finish at time_max
  expect_error(sched$ode_times <- c(0.0, 2, 3), "Last time must be exactly max_time")
  ## Is not sorted:
  expect_error(sched$ode_times <- sched$max_time * c(0, .5, .3, 1), "ode_times must be sorted")
  ## ...and check that none of these caused the times to be set
  expect_false(sched$use_ode_times)
  expect_equal(sched$ode_times, numeric(0))

  ## So, now, manually get the times set up.  For a real challenge, we
  ## should add some more exact hits in here.  Building the expected
  ## times in R is hard enough!
  expected <- rbind(data.frame(species_index=1, start=t1),
                    data.frame(species_index=2, start=t2))
  expected <- expected[order(expected$start, -expected$species_index),]
  expected$end <- c(expected$start[-1], max_t)

  t_ode <- seq(0, sched$max_time, length.out=14)
  idx <- findInterval(t_ode, c(expected$start, max_t), TRUE)
  tmp <- unname(t(apply(expected[c("start", "end")], 1, unlist)))

  expected_ode <- lapply(seq_len(nrow(expected)), function(i)
                         c(tmp[i,1],
                           setdiff(t_ode[idx == i], tmp[i,]),
                           tmp[i,2]))

  ## New schedule because setting and resetting may have changed cohort
  ## order.
  sched <- CohortSchedule(n_species)
  sched$set_times(t1, 1L)
  sched$set_times(t2, 2L)
  sched$max_time <- max_t
  sched$ode_times <- t_ode

  expect_false(sched$use_ode_times)
  expect_identical(sched$ode_times, t_ode)
  sched$use_ode_times <- TRUE
  expect_true(sched$use_ode_times)

  cmp <- drain_schedule(sched)

  expect_equal(sapply(cmp, first), expected$species_index)
  expect_equal(sapply(cmp, second), expected$start)
  expect_equal(sapply(cmp, last), expected$end)

  expect_equal(lapply(cmp, function(x) x[3:(length(x) - 1)]), expected_ode)

  ## check we can clear times:
  sched$clear_ode_times()
  expect_false(sched$use_ode_times)
  expect_equal(sched$ode_times, numeric(0))

  sched$max_time <- Inf
  sched$ode_times <- t_ode
  expect_false(sched$use_ode_times)
  expect_identical(sched$ode_times, t_ode)
  expect_identical(sched$max_time, max(t_ode))
})

test_that("Can expand CohortSchedule", {
  sched <- CohortSchedule(1)
  max_t <- 10
  times1 <- sort(runif(10))
  sched$max_time <- max_t
  sched$set_times(times1, 1)

  expect_equal(sched$n_species, 1)
  expect_identical(sched$max_time, max_t)
  expect_identical(sched$times(1), times1)

  times2 <- sort(runif(20))
  sched3 <- sched$expand(2, times2) # expand by two species:
  expect_equal(sched3$n_species, 3)
  expect_identical(sched3$max_time, max_t)
  expect_identical(sched3$times(1), times1)
  expect_identical(sched3$times(2), times2)
  expect_identical(sched3$times(3), times2)

  # Schedules are independent:
  sched3$max_time <- 2 * max_t
  expect_identical(sched$max_time, max_t)
  expect_identical(sched3$max_time, 2 * max_t)
})
