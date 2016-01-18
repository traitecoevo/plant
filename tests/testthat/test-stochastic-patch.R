context("StochasticPatch")

strategy_types <- get_list_of_strategy_types()

test_that("empty", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi/2,
                          is_resident=TRUE)
    patch <- StochasticPatch(x)(p)

    expect_that(patch, is_a(sprintf("StochasticPatch<%s>",x)))

    expect_that(patch$size, equals(1))
    expect_that(patch$height_max, equals(0.0))
    expect_that(patch$canopy_openness(0), equals(1.0))
    expect_that(patch$ode_state, equals(numeric(0)))
    expect_that(patch$ode_rates, equals(numeric(0)))

    sp <- patch$species
    expect_that(is.list(sp), is_true())
    expect_that(length(sp), equals(1))
    expect_that(sp[[1]], is_a(sprintf("StochasticSpecies<%s>",x)))
    expect_that(sp[[1]]$size, equals(0))
  }
})

test_that("non empty", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi/2,
                          is_resident=TRUE)
    patch <- StochasticPatch(x)(p)
    cmp <- Plant(x)(p$strategies[[1]])

    expect_that(patch$add_seed(0), throws_error("Invalid value"))
    expect_that(patch$add_seed(10), throws_error("out of bounds"))

    expect_that(patch$add_seed(1), is_true())
    expect_that(patch$height_max, is_more_than(0.0))
    expect_that(patch$height_max, equals(cmp$height))

    expect_that(patch$deaths(), equals(0))

    le <- patch$environment$light_environment
    expect_that(range(le$x), equals(c(0.0, cmp$height)))
    expect_that(max(le$y), equals(1.0))
    expect_that(le$y[[1]], is_less_than(1.0))

    if (x == "FF16") {
      expect_that(all(patch$ode_rates > 0.0), is_true())
    } else if (x == "FF16r") {
      expect_that(all(patch$ode_rates[-3] > 0.0), is_true())
      expect_that(patch$ode_rates[[3]], equals(0))
    }
  }
})
