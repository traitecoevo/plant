

test_that("assembly_parameters", {
  expect_is(assembly_parameters(type = "FF16"), "Parameters<FF16,FF16_Env>")

  ## Modify some things:
  p <- assembly_parameters(B_kl2 = 2, type = "FF16")
  p <- assembly_parameters(B_kl2 = 2, a_bio = pi, type = "FF16")
  expect_equal(p$strategy_default$a_bio, pi, type = "FF16")

  expect_error(p2 <- assembly_parameters(list(B_kl2 = 2, a_bio = pi),
    type = "FF16"
  ), "named")
  p2 <- assembly_parameters(pars = list(B_kl2 = 2, a_bio = pi), type = "FF16")
  expect_equal(p2, p)
})

