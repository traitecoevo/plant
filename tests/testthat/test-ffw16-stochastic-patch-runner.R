context("StochasticPatchRunner")

test_that("empty", {
  p <- FFW16_Parameters(strategies=list(FFW16_Strategy()),
                        seed_rain=pi/2,
                        is_resident=TRUE)

  obj <- FFW16_StochasticPatchRunner(p)
  expect_that(obj$time, is_identical_to(0.0))

  sched <- obj$schedule
  expect_that(sched$size, equals(0))
  expect_that(sched$max_time, equals(p$cohort_schedule_max_time))

  ## Now, create a new set of times:
  sched2 <- stochastic_schedule(p)
  expect_that(sched2$size, is_more_than(0))

  obj$schedule <- sched2
  expect_that(obj$schedule$size, equals(sched2$size))

  ## Importantly, this moves time forward to where the first
  ## introduction will be!
  expect_that(obj$time, is_identical_to(sched2$next_event$time_introduction))

  ## I don't think that this is working correctly because the ODE
  ## state is not growing and it should get bigger with each
  ## introduction.
  ##   obj$run_next()
})
