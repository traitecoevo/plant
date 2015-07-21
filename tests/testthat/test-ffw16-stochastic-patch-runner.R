context("StochasticPatchRunner")

test_that("empty", {
  set.seed(1)
  p <- FFW16_Parameters(strategies=list(FFW16_Strategy()),
                        seed_rain=pi/2,
                        is_resident=TRUE,
                        control=fast_control())

  obj <- FFW16_StochasticPatchRunner(p)
  expect_that(obj$time, is_identical_to(0.0))

  sched <- obj$schedule
  expect_that(sched$size, equals(0))
  expect_that(sched$max_time, equals(p$cohort_schedule_max_time))

  ## Now, create a new set of times:
  sched2 <- stochastic_schedule(p)
  expect_that(sched2$size, is_more_than(0))

  ## Does thuis need to happen twice?
  obj$schedule <- sched2
  expect_that(obj$schedule$size, equals(sched2$size))

  ## Importantly, this moves time forward to where the first
  ## introduction will be!
  expect_that(obj$time, is_identical_to(sched2$next_event$time_introduction))

  ## We're empty though....
  expect_that(obj$patch$species[[1]]$size, equals(0))
  expect_that(obj$patch$ode_state, equals(numeric(0)))

  res <- obj$run_next()
  expect_that(res, equals(1L))
  expect_that(obj$time, is_identical_to(sched2$all_times[[1]][[2]]))

  expect_that(length(obj$patch$ode_state), equals(3))
  expect_that(obj$patch$size, equals(1))

  expect_that(obj$complete, is_false())
})

test_that("collect", {
  set.seed(1)
  p <- FFW16_Parameters(strategies=list(FFW16_Strategy()),
                        seed_rain=5/50,
                        patch_area=50,
                        is_resident=TRUE,
                        control=fast_control())
  res <- run_stochastic_collect(p)

  ## This shows that we're probably over-agressively killing plants.
  ## Not sure why, but might be mostly due to the patch area being far
  ## too low.
  if (FALSE) {
    image(attr(res$species, "is_alive")[[1]])
    matplot(res$time, res$species[[1]]["height", , ], type="l",
            lty=1, col="#00000055")
  }
})
