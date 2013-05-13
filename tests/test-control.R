source("helper-tree.R")

context("Control")

ctrl <- new(Control)

obj <- ctrl$parameters

expected <- list(
  cohort_gradient_eps = 1e-6,
  cohort_gradient_richardson = as.numeric(FALSE),
  cohort_gradient_richardson_depth = as.numeric(4L)
  )

keys <- sort(names(expected))
expect_that(sort(names(obj)),
            is_identical_to(sort(names(expected))))
expect_that(obj[keys], is_identical_to(expected[keys]))
