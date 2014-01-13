source("helper-tree.R")

context("Support functions")

test_that("schedule.from.times works", {
  tt <- seq(0, pi, length=100)
  sched <- schedule.from.times(tt)
  expect_that(sched$n_species, equals(1))
  expect_that(sched$size, equals(length(tt) - 1))
  expect_that(sched$times(1), is_identical_to(tt[-length(tt)]))
  expect_that(sched$max_time, is_identical_to(tt[ length(tt)]))

  expect_that(schedule.from.times(rev(tt)), throws_error())
  expect_that(schedule.from.times(c(tt[[1]], tt)), throws_error())
  expect_that(schedule.from.times(1), throws_error())

  n <- 3
  sched <- schedule.from.times(tt, n)
  expect_that(sched$n_species, equals(n))
  expect_that(sched$size, equals(n * (length(tt) - 1)))
  for (i in seq_len(n))
    expect_that(sched$times(i), is_identical_to(tt[-length(tt)]))
  expect_that(sched$max_time, is_identical_to(tt[ length(tt)]))
})

test_that("default.schedule behaves correctly", {
  sched <- default.schedule(100, pi)
  expect_that(sched$n_species, equals(1))
  expect_that(sched$size, equals(100))
  expect_that(sched$times(1), equals(seq(0, pi, length=101)[-101]))
  expect_that(sched$max_time, is_identical_to(pi))
})

test_that("cohort.introduction.times behaves correctly", {
  set.seed(1)
  max.t <- 100 + runif(1)
  min.step.size <- 1e-5
  max.step.size <- 2.0
  times <- cohort.introduction.times(max.t)

  expect_that(first(times), is_identical_to(0))
  expect_that(last(times),  is_identical_to(max.t))

  eps <- .Machine$double.eps
  expect_that(max(diff(times)), is_at_most(max.step.size  + eps))
  expect_that(min(diff(times)), is_at_least(min.step.size - eps))
})
