context("SCM support")

test_that("collect / make_patch", {
  p0 <- scm_base_parameters()
  p0$disturbance_mean_interval <- 30.0
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, FALSE)

  res <- run_scm_collect(p1)

  st_113 <- scm_state(113, res)
  p1_113 <- make_patch(st_113, p1)

  expect_that(p1_113$ode_state, equals(unlist(st_113$species)))
  expect_that(p1_113$time, equals(st_113$time))
  expect_that(p1_113$environment$light_environment$xy,
              equals(unname(st_113$light_env)))
  expect_that(p1_113$canopy_openness(0), is_less_than(0.5))
  expect_that(p1_113$height_max,         is_more_than(10))

  cmp_patch_density <-
    Disturbance(p1$disturbance_mean_interval)$density(res$time)
  expect_that(res$patch_density,
              equals(cmp_patch_density))

  dat <- patch_to_internals(p1_113)
  expect_that(dat, is_a("list"))
  expect_that(length(dat), equals(1))

  dat <- dat[[1]]
  expect_that(dat, is_a("matrix"))
  expect_that(nrow(dat), equals(length(p1_113$species[[1]]$cohorts)))
  n_int <- length(PlantPlus("FF16")(p1$strategies[[1]])$internals)
  expect_that(ncol(dat), equals(n_int + 2L))

  ## NOTE: this currently takes *longer* than the SCM to run due to (I
  ## think) the RcppR6 calling being pretty inefficient in this case.
  ## I should really benchmark it and see why it is so slow, but
  ## possibly we could do this within C++ for a major speedup.
  ints <- scm_to_internals(res)

  n_times <- length(p1$cohort_schedule_times[[1]])
  expect_that(length(ints), equals(1))
  ints <- ints[[1]]
  expect_that(ints, is_a("array"))
  expect_that(length(dim(ints)), equals(3))
  expect_that(dim(ints), equals(c(n_int + 2L, n_times + 1, n_times)))

  v <- rownames(res$species[[1]])
  expect_that(all(v %in% rownames(ints)), is_true())
  expect_that(ints[v, , ], equals(res$species[[1]]))
})

test_that("expand_parameters", {
  p0 <- scm_base_parameters()
  p1 <- expand_parameters(trait_matrix(0.1, "lma"), p0, FALSE)
  ## This will trigger rebuilding the times:
  p1$cohort_schedule_max_time <- 100
  expect_that(p2 <- expand_parameters(trait_matrix(0.2, "lma"), p1, FALSE),
              not(throws_error()))
})
