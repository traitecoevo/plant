context("Support")

test_that("Some support functions", {
  expect_that(fast_control(),        is_a("Control"))
  expect_that(equilibrium_verbose(), is_a("Control"))
  expect_that(equilibrium_quiet(),   is_a("Control"))
  p <- ebt_base_parameters()
  expect_that(p, is_a("Parameters"))
  cmp <- equilibrium_verbose(fast_control())
  cmp$schedule_eps=0.005
  expect_that(p$control, equals(cmp))
})
