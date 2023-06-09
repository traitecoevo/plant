context("Support")

test_that("Some support functions", {
  expect_is(fast_control(), "Control")
  expect_is(equilibrium_verbose(), "Control")
  expect_is(equilibrium_quiet(), "Control")
  base <- scm_base_control()
  expect_is(base, "Control")
  cmp <- equilibrium_verbose(fast_control())
  cmp$schedule_eps <- 0.005
  cmp$equilibrium_eps <- 1e-3
  expect_equal(base, cmp)
})
