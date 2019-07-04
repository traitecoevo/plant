context("Internals")

strategy_types <- get_list_of_strategy_types()

test_that("internals getters and setters", {
  n = 3
  ints = Internals(s_size = n)
  for ( i in 0:(n-1)) {
    ints$set_state(i,10)
    expect_equal(ints$state(i),10)
    ints$set_rate(i,i)
    expect_equal(ints$rate(i),i)
  }
})

test_that("Creation and defaults", {
  internals = Internals(s_size = 0)
  expect_is(internals, "Internals")
  expect_equal(internals$state_size, 0)
  n = 10
  ints = Internals(s_size = n)
  expect_equal(all(is.na(ints$rates)),TRUE)
  expect_identical(ints$states, rep(0.0, n))
})

test_that("Resize", {
  internals = Internals(s_size = 0)
  expect_equal(internals$state_size, 0)
  internals$resize(new_size = 20)
  expect_equal(internals$state_size, 20)
  expect_equal(internals$states, rep(0.0, 20))
  expect_equal(all(is.na(ints$rates)),TRUE)
})
