context("Support")

test_that("Some support functions", {
  expect_is(fast_control(), "Control")
  base <- scm_base_control()
  expect_is(base, "Control")
  cmp <- fast_control()
  cmp$schedule_eps <- 0.005
  expect_equal(base, cmp)
})
