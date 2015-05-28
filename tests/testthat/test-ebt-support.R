context("ebt support")

test_that("collect / make_patch", {
  p0 <- ebt_base_parameters()
  p0$disturbance_mean_interval <- 30.0
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, FALSE)

  res <- run_ebt_collect(p1)

  st_113 <- ebt_state(113, res)
  p1_113 <- make_patch(st_113, p1)

  expect_that(p1_113$ode_state, equals(unlist(st_113$species)))
  expect_that(p1_113$time, equals(st_113$time))
  expect_that(p1_113$canopy_openness(0), is_less_than(0.5))
  expect_that(p1_113$height_max,         is_more_than(10))


  cmp_patch_density <-
    Disturbance(p1$disturbance_mean_interval)$density(res$time)
  expect_that(res$patch_density,
              equals(cmp_patch_density))
})
