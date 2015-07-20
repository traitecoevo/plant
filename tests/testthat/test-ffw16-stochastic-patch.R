context("FFW16_StochasticPatch")

test_that("empty", {
  p <- FFW16_Parameters(strategies=list(FFW16_Strategy()),
                        seed_rain=pi/2,
                        is_resident=TRUE)
  patch <- FFW16_StochasticPatch(p)

  expect_that(patch, is_a("FFW16_StochasticPatch"))

  expect_that(patch$size, equals(1))
  expect_that(patch$height_max, equals(0.0))
  expect_that(patch$canopy_openness(0), equals(1.0))
  expect_that(patch$ode_state, equals(numeric(0)))
  expect_that(patch$ode_rates, equals(numeric(0)))

  sp <- patch$species
  expect_that(is.list(sp), is_true())
  expect_that(length(sp), equals(1))
  expect_that(sp[[1]], is_a("FFW16_StochasticSpecies"))
  expect_that(sp[[1]]$size, equals(0))
})

test_that("non empty", {
  p <- FFW16_Parameters(strategies=list(FFW16_Strategy()),
                        seed_rain=pi/2,
                        is_resident=TRUE)
  patch <- FFW16_StochasticPatch(p)
  cmp <- FFW16_Plant(p$strategies[[1]])

  expect_that(patch$add_seed(0), throws_error("Invalid value"))
  expect_that(patch$add_seed(10), throws_error("out of bounds"))

  patch$add_seed(1)
  expect_that(patch$height_max, is_more_than(0.0))
  expect_that(patch$height_max, equals(cmp$height))

  expect_that(patch$deaths(), equals(0))
})
