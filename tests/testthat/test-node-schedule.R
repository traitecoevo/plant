# Walk the schedule the way a run does: one entry per introduction, whatever
# species it names. There is no second shape for a schedule holding a recording,
# because an introduction is the same thing either way.
drain_schedule <- function(sched) {
  sched$reset()
  cmp <- vector("list", sched$size)
  for (i in seq_len(sched$size)) {
    cmp[[i]] <- sched$next_introduction
    sched$pop()
  }
  cmp
}

drain_column <- function(cmp, field) {
  vapply(cmp, function(x) x[[field]], numeric(1))
}


test_that("NodeScheduleEvent", {
  e <- NodeScheduleEvent(pi, 1)

  expect_identical(e$target_index, 1L)
  expect_identical(e$species_index_raw, 0.0)

  e$target_index <- 2L
  expect_identical(e$target_index, 2L)
  expect_identical(e$species_index_raw, 1.0)

  expect_identical(e$times, pi)
  expect_identical(e$time_introduction, pi)
  expect_identical(e$time_end, pi)
})

test_that("Empty NodeSchedule", {
  n_species <- 2
  sched <- plant:::NodeSchedule(n_species)

  expect_inherits(sched, "NodeSchedule")
  expect_equal(sched$size, 0)
  expect_equal(sched$n_species, n_species)

  expect_equal(sched$remaining, 0)
  expect_equal(sched$max_time, Inf)
  expect_error(sched$next_introduction)
  expect_false(sched$using_ode_steps)
  expect_equal(sched$ode_times, numeric(0))
})

test_that("Corner cases", {
  n_species <- 2
  sched <- plant:::NodeSchedule(n_species)

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
  sched <- plant:::NodeSchedule(n_species)
  sched$set_times(t1, 1L)

  expect_equal(sched$size, length(t1))
  expect_equal(sched$remaining, length(t1))

  species_index <- 1
  expect_equal(sched$times(species_index), t1)
  expect_false(sched$using_ode_steps)
  e <- sched$next_introduction
  expect_identical(e$time, t1[[1]])
  expect_equal(e$species, species_index)

  cmp <- drain_schedule(sched)

  expect_equal(sched$remaining, 0)
  expect_equal(drain_column(cmp, "time"), t1)
  # One interval ends where the next begins, and the last at max_time.
  expect_equal(drain_column(cmp, "time_end"), c(t1[-1], Inf))
  # One species set, so every introduction names it and only it.
  for (e in cmp) {
    expect_equal(e$species, species_index)
  }
})

test_that("Set times (two species)", {
  set.seed(1)
  t1 <- sort(c(0.0, runif(10)))
  t2 <- sort(c(0.0, runif(12)))

  n_species <- 2
  sched <- plant:::NodeSchedule(n_species)
  sched$set_times(t1, 1L)
  sched$set_times(t2, 2L)

  expect_equal(sched$times(1), t1)
  expect_equal(sched$times(2), t2)
  ## Introductions, not species-times: both species start at zero, so that time
  ## is one introduction naming both of them.
  expected_times <- sort(unique(c(t1, t2)))
  expect_equal(sched$size, length(expected_times))
  expect_lt(sched$size, length(t1) + length(t2))

  ## Force a max_time for this run through:
  max_t <- max(c(t1, t2)) + mean(diff(sort(c(t1, t2))))
  sched$max_time <- max_t

  ## Come up with the expected times. Species order within a tied time follows
  ## insertion order now (#628 made the queue's tie-break stable, so two events
  ## of one type at one instant keep the order they were given). Species 1's
  ## times are set first, so at a tie species 1 comes first -- this used to be
  ## the other way round.
  expected <- rbind(data.frame(species_index=1, start=t1),
                    data.frame(species_index=2, start=t2))
  expected <- expected[order(expected$start, expected$species_index),]
  expected$end <- c(expected$start[-1], max_t)

  cmp <- drain_schedule(sched)

  expect_equal(drain_column(cmp, "time"), expected_times)
  expect_equal(drain_column(cmp, "time_end"),
               c(expected_times[-1], max_t))

  ## Species ascending within an introduction, which is the order
  ## Patch::introduced_at rebuilds off the same times -- so the schedule and the
  ## patch agree rather than being each other's reverse.
  for (i in seq_along(cmp)) {
    expect_equal(cmp[[i]]$species,
                 c(if (expected_times[[i]] %in% t1) 1,
                   if (expected_times[[i]] %in% t2) 2))
  }
  expect_equal(cmp[[1]]$species, c(1, 2))

  expect_error(sched$next_introduction, "All introductions completed")
  expect_equal(sched$max_time, max_t)
  sched$reset()
  expect_identical(sched$next_introduction$time, min(c(t1, t2)))
})

test_that("Resetting times replaces them", {
  set.seed(1)
  t1 <- sort(c(0.0, runif(10)))
  t2 <- sort(c(0.0, runif(12)))

  n_species <- 2
  sched <- plant:::NodeSchedule(n_species)
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

  sched <- plant:::NodeSchedule(2)
  sched$set_times(t1, 1)

  last_event <- function(x) {
    x <- x$copy()
    while (x$remaining > 1L) {
      x$pop()
    }
    x$next_introduction
  }

  ## Before setting max_time, the finishing time will be Inf:
  e <- last_event(sched)
  expect_equal(e$time, dplyr::last(t1))
  expect_equal(e$time_end, Inf)

  ## Set max_time to something stupid:
  expect_error(sched$max_time <- 0.5, "max_time must be at least the final")
  expect_error(sched$max_time <- max(t1) - 1e-8, "max_time must be at least the final")

  ## And to something sensible:
  max_t <- max(t1) + 0.1
  sched$max_time <- max_t

  ## Make sure that the last introduction has been modified:
  e <- last_event(sched)
  expect_equal(e$time, dplyr::last(t1))
  expect_equal(e$time_end, max_t)

  ## Now this will fail
  expect_error(sched$set_times(t1 * 2, 1), "Times cannot be greater than max_time")
})

test_that("max_time can be set on an empty schedule", {
  ## An empty schedule constrains max_time not at all, and this is the normal
  ## case: make_node_schedule() and node_schedule_default() both set max_time
  ## before adding any times. The bound is read off the last event, so the
  ## empty case has to be tested explicitly -- reading it unguarded was
  ## undefined behaviour that happened not to crash.
  sched <- plant:::NodeSchedule(1)
  expect_equal(sched$size, 0)
  expect_silent(sched$max_time <- 10)
  expect_equal(sched$max_time, 10)

  ## Still nonnegative-only, and the bound reappears once there are events.
  expect_error(sched$max_time <- -1, "max_time must be nonnegative")
  sched$set_times(c(0, 5), 1)
  expect_error(sched$max_time <- 1, "max_time must be at least the final")
})

test_that("Bulk get/set of times works", {
  n <- 3
  sched <- plant:::NodeSchedule(n)

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

  sched <- plant:::NodeSchedule(n_species)
  sched$set_times(t1, 1L)
  sched$set_times(t2, 2L)
  sched$max_time <- max_t

  none <- numeric(0)
  sched$set_ode_steps(none, none)
  expect_identical(sched$ode_times, none)

  ## Too few values:
  expect_error(sched$set_ode_steps(c(0.0), none), "Need at least two times")
  ## Does not start at 0
  expect_error(sched$set_ode_steps(c(1, 2, 3), none),
               "First time must be exactly zero")
  ## Does not finish at time_max
  expect_error(sched$set_ode_steps(c(0.0, 2, 3), none),
               "Last time must be exactly max_time")
  ## Is not sorted:
  expect_error(sched$set_ode_steps(sched$max_time * c(0, .5, .3, 1), none),
               "ode_times must be sorted")
  ## ...and check that none of these caused the times to be set
  expect_false(sched$using_ode_steps)
  expect_equal(sched$ode_times, numeric(0))

  ## So, now, manually get the times set up.  For a real challenge, we
  ## should add some more exact hits in here.  Building the expected
  ## times in R is hard enough!
  expected <- rbind(data.frame(species_index=1, start=t1),
                    data.frame(species_index=2, start=t2))
  ## Ascending species index at a tie: the queue's tie-break is stable now, so
  ## events keep insertion order (#628).
  expected <- expected[order(expected$start, expected$species_index),]
  expected$end <- c(expected$start[-1], max_t)

  t_ode <- seq(0, sched$max_time, length.out=14)

  ## New schedule because setting and resetting may have changed node
  ## order.
  sched <- plant:::NodeSchedule(n_species)
  sched$set_times(t1, 1L)
  sched$set_times(t2, 2L)
  sched$max_time <- max_t
  # A grid: times with no sizes, so the solver steps to each of them. A schedule
  # holding one uses it -- there is nothing to switch on.
  sched$set_ode_steps(t_ode, numeric(0))
  expect_true(sched$using_ode_steps)
  expect_identical(sched$ode_times, t_ode)
  expect_true(all(is.na(sched$ode_step_sizes)))

  ## Installing a recording does not change what the introductions are: the steps
  ## are not copied into each interval on reset, and an interval's own time list
  ## does not grow to match.
  cmp <- drain_schedule(sched)
  expect_equal(drain_column(cmp, "time"), sort(unique(c(t1, t2))))
  expect_equal(dplyr::last(drain_column(cmp, "time_end")), max_t)

  ## Which recorded steps fall inside which interval is no longer stored, so
  ## there is nothing here that can be stale. It is arithmetic at the point of
  ## use, and that a pinned replay reproduces its own run is checked end to end
  ## in test-scm.R.

  ## check we can clear times:
  sched$clear_ode_steps()
  expect_false(sched$using_ode_steps)
  expect_equal(sched$ode_times, numeric(0))

  sched$max_time <- Inf
  sched$set_ode_steps(t_ode, numeric(0))
  expect_true(sched$using_ode_steps)
  expect_identical(sched$ode_times, t_ode)
  expect_identical(sched$max_time, max(t_ode))
})

test_that("Can expand NodeSchedule", {
  sched <- plant:::NodeSchedule(1)
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
