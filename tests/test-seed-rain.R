source("helper-tree.R")

context("Seed Rain")

x <- c(.1, .2, .3)
rain <- new(SeedRain, length(x))
rain$seed_rain <- x

expect_that(rain$size, equals(length(x)))

rain$first()
expect_that(rain$get(), is_identical_to(x[[1]])) # 1/3
rain$"next"()
expect_that(rain$get(), is_identical_to(x[[2]])) # 2/3
rain$"next"()
expect_that(rain$get(), is_identical_to(x[[3]])) # 3/3
rain$"next"()
expect_that(rain$get(), throws_error())          # 4/3

expect_that(rain$seed_rain <- c(.1, .2, .3, .4), throws_error())

y <- x * 2
rain$seed_rain <- y
rain$first()
expect_that(rain$get(), is_identical_to(y[[1]])) # 1/3

rain2 <- seed_rain(x)
expect_that(rain2$seed_rain, is_identical_to(x))
