context("StochasticPatchRunner")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("empty", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    set.seed(1)
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi/2,
                          is_resident=TRUE,
                          control=fast_control())
    
    if(grepl("K93", x))
      p$k_I <- 1e-3

    obj <- StochasticPatchRunner(x, e)(p)
    expect_identical(obj$time, 0.0)

    sched <- obj$schedule
    expect_equal(sched$size, 0)
    expect_equal(sched$max_time, p$cohort_schedule_max_time)

    ## Now, create a new set of times:
    sched2 <- stochastic_schedule(p)
    expect_gt(sched2$size, 0)

    ## Does thuis need to happen twice?
    obj$schedule <- sched2
    expect_equal(obj$schedule$size, sched2$size)

    ## Importantly, this moves time forward to where the first
    ## introduction will be!
    expect_identical(obj$time, sched2$next_event$time_introduction)

    ## We're empty though....
    expect_equal(obj$patch$species[[1]]$size, 0)
    expect_equal(obj$patch$ode_state, numeric(0))

    res <- obj$run_next()
    expect_equal(res, 1L)
    expect_identical(obj$time, sched2$all_times[[1]][[2]])

    ode_size <- Individual(x, e)(strategy_types[[x]]())$ode_size
    expect_equal(length(obj$patch$ode_state), ode_size)
    expect_equal(obj$patch$size, 1)

    expect_false(obj$complete)
  }
})

test_that("collect", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    set.seed(1)
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          seed_rain=5/50,
                          patch_area=50,
                          is_resident=TRUE,
                          control=fast_control())
    if(grepl("K93", x))
      p$k_I <- 1e-3
    
    expect_silent(res <- run_stochastic_collect(p))
    ## TODO: more tests on collect output

    ## This shows that we're probably over-aggressively killing plants.
    ## Not sure why, but might be mostly due to the patch area being far
    ## too low.
    if (FALSE) {
      image(attr(res$species, "is_alive")[[1]])
      matplot(res$time, res$species[[1]]["height", , ], type="l",
              lty=1, col="#00000055")
    }
  }
})
