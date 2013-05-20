source("helper-tree.R")

context("Seed Rain")

x <- c(.1, .2, .3)
rain <- new(SeedRain, x)
rain$first()
expect_that(rain$get(), is_identical_to(x[[1]])) # 1/3
rain$"next"()
expect_that(rain$get(), is_identical_to(x[[2]])) # 2/3
rain$"next"()
expect_that(rain$get(), is_identical_to(x[[3]])) # 3/3
rain$"next"()
expect_that(rain$get(), throws_error())          # 4/3

expect_that(rain$set(c(.1, .2, .3, .4)), throws_error())

y <- x * 2
rain$set(y)
rain$first()
expect_that(rain$get(), is_identical_to(y[[1]])) # 1/3
