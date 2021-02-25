context("Internals")

test_that("internals getters and setters", {
  n = 3
  a_n = 2
  ints = Internals(s_size = n, a_size = 2)
  for ( i in 0:(n-1)) {
    ints$set_state(i,10)
    expect_equal(ints$state(i), 10)
    ints$set_rate(i, i)
    expect_equal(ints$rate(i), i)
  }

  expect_equal(length(ints$auxs), a_n)
  expect_equal(length(ints$states), n)
  expect_equal(length(ints$rates), length(ints$states))
})

test_that("Creation and defaults", {
  internals = Internals(s_size = 0, a_size = 0)
  expect_is(internals, "Internals")
  expect_equal(internals$state_size, 0)
  expect_equal(internals$aux_size, 0)
  n = 10
  ints = Internals(s_size = n, a_size = n)
  expect_equal(all(is.na(ints$rates)),TRUE)
  expect_identical(ints$states, rep(0.0, n))
})

test_that("Resize", {
  internals = Internals(s_size = 0, a_size = 0)
  expect_equal(internals$state_size, 0)
  expect_equal(internals$aux_size, 0)
  internals$resize(new_size = 20, new_aux_size = 10)
  expect_equal(internals$state_size, 20)
  expect_equal(internals$states, rep(0.0, 20))
  expect_equal(internals$auxs, rep(0.0, 10))
  expect_equal(all(is.na(internals$rates)),TRUE)
})
