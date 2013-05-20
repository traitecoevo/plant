source("helper-tree.R")

context("Seed Rain")

rain <- new(SeedRain, 3)
rain$first()
expect_that(rain$get(), is_identical_to(0.0)) # 1/3
rain$"next"()
expect_that(rain$get(), is_identical_to(0.0)) # 2/3
rain$"next"()
expect_that(rain$get(), is_identical_to(0.0)) # 3/3
rain$"next"()
expect_that(rain$get(), throws_error())       # 4/3

rain$set(c(.1, .2, .3))
expect_that(rain$set(c(.1, .2, .3, .4)), throws_error())

rain$first()
expect_that(rain$get(), is_identical_to(0.1)) # 1/3
rain$"next"()
expect_that(rain$get(), is_identical_to(0.2)) # 2/3
rain$"next"()
expect_that(rain$get(), is_identical_to(0.3)) # 3/3
