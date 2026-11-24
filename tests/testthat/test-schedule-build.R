
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("refine_schedule on no residents runs trivially", {
  for (x in names(strategy_types)) {
    p <- testing_parameters(x)
    expect_silent(solver <- run(p, refine_schedule = TRUE)$solver)
    expect_equal(length(solver$parameters$strategies), 0)
  }
})

test_that("bisection leaves an interval whose ends are adjacent doubles whole", {
  ## On [1, 2) adjacent doubles are .Machine$double.eps apart, and the
  ## midpoint of two of them rounds onto one end.
  a <- 1.25
  b <- a + .Machine$double.eps
  expect_true(b > a)
  expect_true((0.5 * (a + b)) %in% c(a, b))

  times <- c(0, a, b, 2)
  res <- bisect_flagged_intervals(times, c(FALSE, TRUE, TRUE, TRUE))
  expect_identical(res$times, sort(c(times, 0.5 * a, 0.5 * (b + 2))))
  expect_identical(res$unsplit, b)
  expect_identical(anyDuplicated(res$times), 0L)

  ## Every other interval is bisected at its midpoint, as it always was.
  times <- c(0, 0.5, 1.5, 4)
  res <- bisect_flagged_intervals(times, c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(res$times, c(0, 0.25, 0.5, 1.5, 2.75, 4))
  expect_identical(res$unsplit, numeric(0))
})

test_that("Schedule building", {
  for (x in c("FF16")) {
    ## This is a really dumb test but it should act as a regression test
    ## at least. Refinement now happens entirely in C++ via
    ## run(refine_schedule = TRUE) -> DeterministicPatchRunner::refine_schedule().
    p <- testing_parameters(x)
    p$strategies <- list(strategy_types[[x]]())
    p$strategies[[1]]$birth_rate_y <- 0.1

    env <- Environment(x)
    ctrl <- Control()

    solver <- run(p, env, ctrl, refine_schedule = TRUE)$solver
    p2 <- solver$parameters

    expect_equal(length(p2$node_schedule_times_default), 141)
    expect_equal(length(p2$node_schedule_times[[1]]), 148)

    ## Refinement leaves the Parameters self-describing (schedule + ode times).
    expect_equal(length(p2$ode_times), length(solver$ode_times))
    expect_true(length(solver$patch$offspring_production) == 1L)
  }
})
