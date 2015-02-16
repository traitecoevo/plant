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

  expect_that(e$species_index, is_identical_to(1L))
  expect_that(e$species_index_raw, is_identical_to(0.0))

  e$species_index <- 2L
  expect_that(e$species_index, is_identical_to(2L))
  expect_that(e$species_index_raw, is_identical_to(1.0))

  expect_that(e$times,             is_identical_to(pi))
  expect_that(e$time_introduction, is_identical_to(pi))
  expect_that(e$time_end,          is_identical_to(pi))
})

test_that("Empty CohortSchedule", {
  n_species <- 2
  sched <- CohortSchedule(n_species)

  expect_that(sched, is_a("CohortSchedule"))
  expect_that(sched$size,          equals(0))
  expect_that(sched$n_species,     equals(n_species))

  expect_that(sched$remaining,     equals(0))
  expect_that(sched$max_time,      equals(Inf))
  expect_that(sched$next_event,    throws_error())
  expect_that(sched$use_ode_times, is_false())
  expect_that(sched$ode_times,     equals(numeric(0)))
})

test_that("Corner cases", {
  n_species <- 2
  sched <- CohortSchedule(n_species)

  set.seed(1)
  t1 <- c(0.0, runif(10))
  t2 <- c(0.0, runif(12))

  expect_that(sched$set_times(t1, 1L),
              throws_error("Times must be sorted"))

  t1 <- sort(t1)
  t2 <- sort(t2)

  expect_that(sched$set_times(t1, -1),
              throws_error("Invalid value for index"))
  expect_that(sched$set_times(t1, 0),
              throws_error("Invalid value for index"))
  expect_that(sched$set_times(t1, n_species + 1L),
              throws_error("Index 3 out of bounds"))
})

test_that("Set times (one species)", {
  set.seed(1)
  t1 <- sort(c(0.0, runif(10)))
  t2 <- sort(c(0.0, runif(12)))

  n_species <- 2
  sched <- CohortSchedule(n_species)
  sched$set_times(t1, 1L)

  expect_that(sched$size,      equals(length(t1)))
  expect_that(sched$remaining, equals(length(t1)))

  species_index <- 1
  expect_that(sched$times(species_index), equals(t1))
  expect_that(sched$use_ode_times, is_false())
  e <- sched$next_event
  expect_that(e$time_introduction,   is_identical_to(t1[[1]]))
  expect_that(e$species_index,       equals(species_index))

  cmp <- drain_schedule(sched)

  expect_that(sched$remaining, equals(0))
  n <- length(t1)
  expect_that(cmp[,1], equals(rep(1, n)))
  expect_that(cmp[,2], equals(t1))
  expect_that(cmp[,3], equals(t1))
  expect_that(cmp[,4], equals(c(t1[-1], Inf)))
  expect_that(cmp[,5], equals(c(t1[-1], Inf)))
})

test_that("Set times (two species)", {
  set.seed(1)
  t1 <- sort(c(0.0, runif(10)))
  t2 <- sort(c(0.0, runif(12)))

  n_species <- 2
  sched <- CohortSchedule(n_species)
  sched$set_times(t1, 1L)
  sched$set_times(t2, 2L)

  expect_that(sched$times(1), equals(t1))
  expect_that(sched$times(2), equals(t2))
  expect_that(sched$size, equals(length(t1) + length(t2)))

  ## Force a max_time for this run through:
  max_t <- max(c(t1, t2)) + mean(diff(sort(c(t1, t2))))
  sched$max_time <- max_t

  ## Come up with the expected times (2nd argument ensures stable sort)
  expected <- rbind(data.frame(species_index=1, start=t1),
                    data.frame(species_index=2, start=t2))
  expected <- expected[order(expected$start, -expected$species_index),]
  expected$end <- c(expected$start[-1], max_t)

  cmp <- drain_schedule(sched)

  expect_that(cmp[,1], equals(expected$species_index))
  expect_that(cmp[,2], equals(expected$start))
  expect_that(cmp[,3], equals(expected$start))
  expect_that(cmp[,4], equals(expected$end))
  expect_that(cmp[,5], equals(expected$end))

  expect_that(sched$next_event, throws_error("All events completed"))
  expect_that(sched$max_time, equals(max_t))
  sched$reset()
  expect_that(sched$next_event$time_introduction,
              is_identical_to(min(c(t1, t2))))
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
  expect_that(sched$times(1), equals(t1_new))
  sched$set_times(t1, 1)
  expect_that(sched$times(1), equals(t1))
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
  expect_that(e$time_introduction, equals(last(t1)))
  expect_that(e$time_end,          equals(Inf))

  ## Set max_time to something stupid:
  expect_that(sched$max_time <- 0.5,
              throws_error("max_time must be at least the final"))
  expect_that(sched$max_time <- max(t1) - 1e-8,
              throws_error("max_time must be at least the final"))

  ## And to something sensible:
  max_t <- max(t1) + 0.1
  sched$max_time <- max_t

  ## Make sure that the last event has been modified:
  e <- last_event(sched)
  expect_that(e$time_introduction, equals(last(t1)))
  expect_that(e$time_end,          equals(max_t))

  ## Now this will fail
  expect_that(sched$set_times(t1 * 2, 1),
              throws_error("Times cannot be greater than max_time"))
})

test_that("Bulk get/set of times works", {
  n <- 3
  sched <- CohortSchedule(n)

  set.seed(1)
  t_new <- lapply(seq_len(n), function(...) sort(runif(rpois(1, 10))))
  sched$all_times <- t_new
  expect_that(sched$all_times, is_identical_to(t_new))

  expect_that(sched$all_times <- t_new[1],
              throws_error("Incorrect length input"))
  expect_that(sched$all_times <- t_new[[1]],
              throws_error("Incorrect length input"))
  expect_that(sched$all_times <- t_new[c(1, seq_len(n))],
              throws_error("Incorrect length input"))
  ## Unfortunately, this does not throw.
  ## expect_that(sched$all_times <- seq_len(n), throws_error())
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
  expect_that(sched$ode_times, is_identical_to(numeric(0)))

  ## Too few values:
  expect_that(sched$ode_times <- c(0.0),
              throws_error("Need at least two times"))
  ## Does not start at 0
  expect_that(sched$ode_times <- c(1, 2, 3),
              throws_error("First time must be exactly zero"))
  ## Does not finish at time_max
  expect_that(sched$ode_times <- c(0.0, 2, 3),
              throws_error("Last time must be exactly max_time"))
  ## Is not sorted:
  expect_that(sched$ode_times <- sched$max_time * c(0, .5, .3, 1),
              throws_error("ode_times must be sorted"))
  ## ...and check that none of these caused the times to be set
  expect_that(sched$use_ode_times, is_false())
  expect_that(sched$ode_times,     equals(numeric(0)))

  ## So, now, manually get the times set up.  For a real challenge, we
  ## should add some more exact hits in here.  Building the expected
  ## times in R is hard enough!
  expected <- rbind(data.frame(species_index=1, start=t1),
                    data.frame(species_index=2, start=t2))
  expected <- expected[order(expected$start, -expected$species_index),]
  expected$end <- c(expected$start[-1], max_t)

  t_ode <- seq(0, sched$max_time, length=14)
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

  expect_that(sched$use_ode_times, is_false())
  expect_that(sched$ode_times,     is_identical_to(t_ode))
  sched$use_ode_times <- TRUE
  expect_that(sched$use_ode_times, is_true())

  cmp <- drain_schedule(sched)

  expect_that(sapply(cmp, first),  equals(expected$species_index))
  expect_that(sapply(cmp, second), equals(expected$start))
  expect_that(sapply(cmp, last),   equals(expected$end))

  expect_that(lapply(cmp, function(x) x[3:(length(x) - 1)]),
              equals(expected_ode))

  ## check we can clear times:
  sched$clear_ode_times()
  expect_that(sched$use_ode_times, is_false())
  expect_that(sched$ode_times,   equals(numeric(0)))

  sched$max_time <- Inf
  sched$ode_times <- t_ode
  expect_that(sched$use_ode_times, is_false())
  expect_that(sched$ode_times,   is_identical_to(t_ode))
  expect_that(sched$max_time,    is_identical_to(max(t_ode)))
})

test_that("Can expand CohortSchedule", {
  sched <- CohortSchedule(1)
  max_t <- 10
  times1 <- sort(runif(10))
  sched$max_time <- max_t
  sched$set_times(times1, 1)

  expect_that(sched$n_species, equals(1))
  expect_that(sched$max_time, is_identical_to(max_t))
  expect_that(sched$times(1), is_identical_to(times1))

  times2 <- sort(runif(20))
  sched3 <- sched$expand(2, times2) # expand by two species:
  expect_that(sched3$n_species, equals(3))
  expect_that(sched3$max_time, is_identical_to(max_t))
  expect_that(sched3$times(1), is_identical_to(times1))
  expect_that(sched3$times(2), is_identical_to(times2))
  expect_that(sched3$times(3), is_identical_to(times2))

  # Schedules are independent:
  sched3$max_time <- 2 * max_t
  expect_that(sched$max_time,  is_identical_to(max_t))
  expect_that(sched3$max_time, is_identical_to(2 * max_t))
})
