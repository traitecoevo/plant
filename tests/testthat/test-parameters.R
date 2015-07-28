
## TODO: Remove ["FFW16"] to test this with all types. But first
## requires issue #162 to be resolved
strategy_types <- get_list_of_strategy_types()["FFW16"]
hyperpar_functions <-get_list_of_hyperpar_functions()

for (x in names(strategy_types)) {
  context(sprintf("Parameters-%s",x))
  test_that("Creation & defaults", {
    s <- strategy_types[[x]]()
    p <- Parameters(x)()
    expect_that(p, is_a(sprintf("Parameters<%s>",x)))

    expect_that(length(p$strategies), equals(0))
    expect_that(length(p$is_resident), equals(0))
    expect_that(length(p$seed_rain), equals(0))

    expected <- list(c_ext=0.5,
                     n_patches=1,    # NOTE: Different to tree 0.1
                     patch_area=1.0, # NOTE: Different to tree 0.1
                     disturbance_mean_interval=30.0,
                     hyperpar=hyperpar_functions[[x]])

    expect_that(p[names(expected)], is_identical_to(expected))
    expect_that(p$strategy_default, equals(s))

    expect_that(p$cohort_schedule_max_time,
                is_identical_to(cohort_schedule_max_time_default(p)))
    expect_that(p$cohort_schedule_times_default,
                is_identical_to(cohort_schedule_times_default(p$cohort_schedule_max_time)))
    expect_that(p$cohort_schedule_times, is_identical_to(list()))
  })

  test_that("Nontrivial creation", {
    s <- strategy_types[[x]]()
    p <- Parameters(x)(strategies=list(s))
    expect_that(p$seed_rain, equals(1.0))
    expect_that(p$is_resident, is_true())

    expect_that(Parameters(x)(seed_rain=pi),
                throws_error("Incorrect length seed_rain"))
    expect_that(Parameters(x)(is_resident=FALSE),
                throws_error("Incorrect length is_resident"))

    expect_that(Parameters(x)(strategies=list(strategy_types[[x]](),
                                   strategy_types[[x]]()),
                                 seed_rain=pi),
                throws_error("Incorrect length seed_rain"))
    expect_that(Parameters(x)(strategies=list(strategy_types[[x]](),
                                   strategy_types[[x]]()),
                                 is_resident=TRUE),
                throws_error("Incorrect length is_resident"))

    p <- Parameters(x)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi,
                          is_resident=TRUE)

    expect_that(p$cohort_schedule_max_time,
                is_identical_to(cohort_schedule_max_time_default(p)))
    expect_that(p$cohort_schedule_times_default,
                is_identical_to(cohort_schedule_times_default(p$cohort_schedule_max_time)))
    expect_that(p$cohort_schedule_times,
                is_identical_to(list(p$cohort_schedule_times_default)))

    ## Now, with some of these set:
    p <- Parameters(x)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi,
                          is_resident=TRUE,
                          disturbance_mean_interval=2)

    expect_that(p$cohort_schedule_max_time,
                is_less_than(10))
    expect_that(p$cohort_schedule_times_default,
                is_identical_to(cohort_schedule_times_default(p$cohort_schedule_max_time)))
    expect_that(p$cohort_schedule_times,
                is_identical_to(list(p$cohort_schedule_times_default)))
  })

  test_that("Parameters overwrites Strategy control", {
    ctrl <- ctrl_s <- ctrl_p <- Control()
    ## set these just as markers:
    ctrl_s$schedule_eps <- 1
    ctrl_p$schedule_eps <- 2

    s <- strategy_types[[x]](control=ctrl_s)
    expect_that(s$control, is_identical_to(ctrl_s))
    expect_that(s$control, not(is_identical_to(ctrl_p)))

    p <- Parameters(x)(control=ctrl_p)
    expect_that(p$control, not(is_identical_to(ctrl_s)))
    expect_that(p$control, is_identical_to(ctrl_p))

    p$strategies <- list(s)
    p$seed_rain <- 1
    p$is_resident <- TRUE
    ## Pass though to force validation:
    tmp <- Patch(x)(p)$parameters
    expect_that(tmp$control, is_identical_to(ctrl_p))
    expect_that(tmp$strategies[[1]]$control, is_identical_to(ctrl_p))

    ## In one shot:
    p2 <- Parameters(x)(control=ctrl_p, strategies=list(s),
                           seed_rain=1, is_resident=TRUE)
    expect_that(p2$control, is_identical_to(ctrl_p))
    expect_that(p2$strategies[[1]]$control, is_identical_to(ctrl_p))
  })

  test_that("Generate cohort schedule", {
    p <- Parameters(x)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi/2, is_resident=TRUE)
    sched <- make_cohort_schedule(p)

    expect_that(sched$n_species, equals(1))
    expect_that(sched$max_time, equals(p$cohort_schedule_max_time))
    expect_that(sched$times(1), equals(p$cohort_schedule_times_default))
    expect_that(sched$all_times, equals(p$cohort_schedule_times))
  })

  test_that("Store hyperparams", {
    p <- Parameters(x)(hyperpar=hyperpar_functions[[x]])
    expect_that(p$hyperpar, is_identical_to(hyperpar_functions[[x]]))
    tmp <- Patch(x)(p)$parameters
    expect_that(tmp$hyperpar, is_identical_to(hyperpar_functions[[x]]))
  })

  test_that("ebt_base_parameters", {
    p <- ebt_base_parameters()
    expect_that(p$hyperpar, equals(hyperpar_functions[[x]]))
  })

  test_that("Disturbance interval", {
    p <- ebt_base_parameters()
    expect_that(p$disturbance_mean_interval, equals(30.0))
    expect_that(p$cohort_schedule_max_time,
                is_identical_to(cohort_schedule_max_time_default(p)))
    p$strategies <- list(strategy_types[[x]]())

    p$disturbance_mean_interval <- 10.0
    ## This is going to force us back through the validator
    p2 <- validate(p)
    expect_that(p2$cohort_schedule_max_time,
                is_identical_to(cohort_schedule_max_time_default(p2)))
    expect_that(p2$cohort_schedule_max_time,
                is_less_than(p$cohort_schedule_max_time))
    expect_that(last(p2$cohort_schedule_times_default),
                is_less_than(p2$cohort_schedule_max_time))
    expect_that(last(p2$cohort_schedule_times_default),
                is_less_than(p2$cohort_schedule_max_time))
    expect_that(p2$cohort_schedule_times,
                equals(list(p2$cohort_schedule_times_default)))

    ## We will blow away any data that is stored in p$cohort_schedule*
    expect_that(p$cohort_schedule_max_time, is_more_than(100))
    p$cohort_schedule_max_time <- 1
    p$cohort_schedule_times_default <- 1:10
    p$cohort_schedule_time <- list(1:11)
    expect_that(validate(p), equals(p2))
  })
}