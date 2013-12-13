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
})

test_that("default.schedule behaves correctly", {
  sched <- default.schedule(100, pi)
  expect_that(sched$n_species, equals(1))
  expect_that(sched$size, equals(100))
  expect_that(sched$times(1), equals(seq(0, pi, length=101)[-101]))
  expect_that(sched$max_time, is_identical_to(pi))
})
