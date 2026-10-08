
test_that("internals getters and setters", {
  n = 3
  a_n = 2
  ints = plant:::Internals(s_size = n, a_size = 2)
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
  internals = plant:::Internals(s_size = 0, a_size = 0)
  expect_inherits(internals, "Internals")
  expect_equal(internals$state_size, 0)
  expect_equal(internals$aux_size, 0)
  n = 10
  ints = plant:::Internals(s_size = n, a_size = n)
  expect_equal(all(is.na(ints$rates)),TRUE)
  expect_identical(ints$states, rep(0.0, n))
})

# Per-layer water uptake is copied out to R with the rest of the store, but had
# no accessor, so the whole water channel of the rate vector was unreadable.
test_that("consumption rates are readable", {
  ints <- plant:::Internals(s_size = 3, a_size = 2)
  expect_equal(ints$resource_size, 0)
  expect_equal(length(ints$consumption_rates), 0)
})

test_that("consumption rates carry one water uptake per soil layer", {
  theta <- rep(0.25, 5)
  ind <- TF24_Individual(TF24_Strategy())
  ind$set_state("height", 5)
  env <- TF24_Environment()
  env$extrinsic_drivers_set_constant("PPFD", 1800)
  env$set_soil_number_of_depths(length(theta))
  env$set_soil_water_state(theta)
  ind$compute_rates(env)

  ints <- ind$internals
  expect_equal(ints$resource_size, length(theta))
  expect_equal(length(ints$consumption_rates), length(theta))
  expect_true(all(is.finite(ints$consumption_rates)))
  # The indexed read and the vector are the same store, 0- and 1-based.
  for (i in seq_along(theta)) {
    expect_equal(ints$consumption_rate(i - 1), ints$consumption_rates[[i]])
  }
  # A lit, well-watered plant transpires, so uptake is drawn from somewhere.
  expect_gt(sum(ints$consumption_rates), 0)
})

test_that("Resize", {
  internals = plant:::Internals(s_size = 0, a_size = 0)
  expect_equal(internals$state_size, 0)
  expect_equal(internals$aux_size, 0)
  internals$resize(new_size = 20, new_aux_size = 10)
  expect_equal(internals$state_size, 20)
  expect_equal(internals$states, rep(0.0, 20))
  expect_equal(internals$auxs, rep(0.0, 10))
  expect_equal(all(is.na(internals$rates)),TRUE)
})
