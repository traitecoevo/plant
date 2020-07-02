context("Support")

test_that("Some support functions", {
  expect_is(fast_control(), "Control")
  expect_is(equilibrium_verbose(), "Control")
  expect_is(equilibrium_quiet(), "Control")
  p <- scm_base_parameters()
  expect_is(p, "Parameters<FF16,FF16_Env>")
  cmp <- equilibrium_verbose(fast_control())
  cmp$schedule_eps <- 0.005
  cmp$equilibrium_eps <- 1e-3
  expect_equal(p$control, cmp)
})

test_that("assembly_parameters", {
  expect_is(assembly_parameters(), "Parameters<FF16,FF16_Env>")

  ## Modify some things:
  p <- assembly_parameters(B_kl2=2)
  p <- assembly_parameters(B_kl2=2, a_bio=pi)
  expect_equal(p$strategy_default$a_bio, pi)

  expect_error(p2 <- assembly_parameters(list(B_kl2=2, a_bio=pi)), "named")
  p2 <- assembly_parameters(pars=list(B_kl2=2, a_bio=pi))
  expect_equal(p2, p)
})
