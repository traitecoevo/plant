source("helper-tree.R")

context("Plant utilities")

test_that("Grow plant to given height", {
  # Stopping height of 10:
  height.10 <- function(plant) {
    plant$height - 10
  }

  env <- test.environment(10)
  p <- new(Plant, new(Strategy))
  h0 <- p$height
  p2 <- grow.plant.to.size(p, env, height.10)
  expect_that(p2$height, equals(10, tolerance=1e-5))
  # Original plant unchanged:
  expect_that(p$height, is_identical_to(h0))
})

test_that("Grow plant to given leaf area", {
  area.1 <- function(plant) {
    plant$leaf_area - 1
  }

  env <- test.environment(10)
  p <- new(Plant, new(Strategy))
  p2 <- grow.plant.to.size(p, env, area.1)
  expect_that(p2$leaf_area, equals(1, tolerance=1e-5))
})

test_that("Corner cases", {
  # Starting value smaller than target (so can't grow)
  height.0 <- function(plant) {
    plant$height - 0
  }

  env <- test.environment(10)
  p <- new(Plant, new(Strategy))

  expect_that(grow.plant.to.size(p, env, height.0), throws_error())

  # Finishing value far too big:
  height.Inf <- function(plant) {
    plant$height - Inf
  }
  expect_that(grow.plant.to.size(p, env, height.Inf, 10), throws_error())
})
