context("add_strategies API (#410)")

test_that("add_strategies matches the deprecated expand_parameters", {
  p0 <- scm_base_parameters("FF16")
  tr <- trait_matrix(c(0.08, 0.26), "lma")

  new <- add_strategies(p0, tr, birth_rate = list(20, 20))
  old <- suppressWarnings(
    expand_parameters(tr, p0, birth_rate_list = list(20, 20)))

  expect_identical(new$strategies, old$strategies)
  expect_identical(new$node_schedule_times, old$node_schedule_times)
})

test_that("add_mutant matches the deprecated mutant_parameters", {
  p0 <- scm_base_parameters("FF16")
  p1 <- add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 20)
  tr <- trait_matrix(0.2, "lma")

  new <- add_mutant(p1, tr, birth_rate = 20)
  old <- suppressWarnings(mutant_parameters(tr, p1, birth_rate_list = 20))

  expect_identical(new$strategies, old$strategies)
  expect_equal(length(new$strategies), 1L)   # mutant replaces residents
})

test_that("generate_strategy matches the deprecated strategy_list", {
  p0 <- scm_base_parameters("FF16")
  tr <- trait_matrix(0.08, "lma")
  new <- generate_strategy(p0, tr, birth_rate = 1)
  old <- suppressWarnings(strategy_list(tr, p0, birth_rate_list = 1))
  expect_identical(new, old)
})

test_that("deprecated shims warn", {
  p0 <- scm_base_parameters("FF16")
  tr <- trait_matrix(0.08, "lma")
  expect_warning(expand_parameters(tr, p0, birth_rate_list = 1), "deprecated")
  expect_warning(strategy_list(tr, p0, birth_rate_list = 1), "deprecated")
  expect_warning(mutant_parameters(tr, p0, birth_rate_list = 1), "deprecated")
})

test_that("add_strategies is pipe-friendly", {
  p <- scm_base_parameters("FF16") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  expect_equal(p$strategies[[1]]$pars$lma, 0.0825)
  out <- run_scm(p, Environment("FF16"), Control())
  expect_false(is.null(out))
})
