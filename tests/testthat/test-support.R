context("Support")

test_that("Some support functions", {
  expect_that(fast_control(),        is_a("Control"))
  expect_that(equilibrium_verbose(), is_a("Control"))
  expect_that(equilibrium_quiet(),   is_a("Control"))
  p <- scm_base_parameters()
  expect_that(p, is_a("Parameters<FF16>"))
  cmp <- equilibrium_verbose(fast_control())
  cmp$schedule_eps <- 0.005
  cmp$equilibrium_eps <- 1e-3
  expect_that(p$control, equals(cmp))
})

test_that("assembly_parameters", {
  expect_that(assembly_parameters(), is_a("Parameters<FF16>"))

  ## Modify some things:
  p <- assembly_parameters(B_kl2=2)
  expect_that(environment(p$hyperpar)$B_kl2, equals(2))
  p <- assembly_parameters(B_kl2=2, a_bio=pi)
  expect_that(environment(p$hyperpar)$B_kl2, equals(2))
  expect_that(p$strategy_default$a_bio, equals(pi))

  expect_that(p2 <- assembly_parameters(list(B_kl2=2, a_bio=pi)),
              throws_error("named"))
  p2 <- assembly_parameters(pars=list(B_kl2=2, a_bio=pi))
  expect_that(p2, equals(p))
})
