context("SCM support")


test_that("collect / make_patch", {
  p0 <- scm_base_parameters()
  hyperpar <- make_FF16_hyperpar()
  p0$disturbance_mean_interval <- 30.0
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, hyperpar, FALSE)

  res <- run_scm_collect(p1)

  st_113 <- scm_state(113, res)
  p1_113 <- make_patch(st_113, p1)

  expect_equal(p1_113$ode_state, unlist(st_113$species))
  expect_equal(p1_113$time, st_113$time)
  expect_equal(p1_113$environment$canopy$canopy_interpolator$xy, unname(st_113$env))
  expect_lt(exp(-p0$k_I * p1_113$compute_competition(0)), 0.5)
  expect_gt(p1_113$height_max, 10)

  cmp_patch_density <- Disturbance(p1$disturbance_mean_interval)$density(res$time)
  expect_equal(res$patch_density, cmp_patch_density)

  dat <- patch_to_internals(p1_113)
  expect_is(dat, "list")
  expect_equal(length(dat), 1)

  dat <- dat[[1]]
  expect_is(dat, "matrix")
  expect_equal(nrow(dat), length(p1_113$species[[1]]$cohorts))
  # once for rates, once for states
  n_int <- (Individual("FF16","FF16_Env")(p1$strategies[[1]])$ode_size * 2) + 
    Individual("FF16","FF16_Env")(p1$strategies[[1]])$aux_size
  cat(ncol(dat), n_int + 2L)
  expect_equal(ncol(dat), n_int + 2L)

  # NOTE: this currently takes *longer* than the SCM to run due to (I
  # think) the RcppR6 calling being pretty inefficient in this case.
  # I should really benchmark it and see why it is so slow, but
  # possibly we could do this within C++ for a major speedup.
  ints <- scm_to_internals(res)

  n_times <- length(p1$cohort_schedule_times[[1]])
  expect_equal(length(ints), 1)
  ints <- ints[[1]]
  expect_is(ints, "array")
  expect_equal(length(dim(ints)), 3)
  expect_equal(dim(ints), c(n_int + 2L, n_times + 1, n_times))

  v <- rownames(res$species[[1]])
  expect_true(all(v %in% rownames(ints)))
  expect_equal(ints[v, , ], res$species[[1]])

})

test_that("expand_parameters", {
  hyperpar <- make_FF16_hyperpar()
  p0 <- scm_base_parameters()
  p1 <- expand_parameters(trait_matrix(0.1, "lma"), p0, mutant=FALSE)
  ## This will trigger rebuilding the times:
  p1$cohort_schedule_max_time <- 100
  expect_silent(p2 <- expand_parameters(trait_matrix(0.2, "lma"), p1, mutant=FALSE))
})
