context("FFW16_StochasticSpecies")

test_that("empty", {
  env <- test_environment(3, seed_rain=1.0)
  s <- FFW16_Strategy()
  sp <- FFW16_StochasticSpecies(s)

  expect_that(sp, is_a("FFW16_StochasticSpecies"))
  expect_that(sp$size, equals(0))

  seed <- sp$seed
  expect_that(seed, is_a("Plant"))
  expect_that(seed, is_a("Plant<FFW16>"))

  expect_that(sp$heights, equals(numeric(0)))
  expect_that(sp$height_max, equals(0.0))
  expect_that(sp$species, is_identical_to(NULL))
  expect_that(sp$plants, equals(list()))
  expect_that(sp$area_leaf_above(0), equals(0.0))
  expect_that(sp$ode_size, equals(0))
  expect_that(sp$ode_state, is_identical_to(numeric(0)))
  expect_that(sp$ode_rates, is_identical_to(numeric(0)))
})

test_that("Single individual", {
  env <- test_environment(3, seed_rain=1.0)
  s <- FFW16_Strategy()
  sp <- FFW16_StochasticSpecies(s)
  p <- FFW16_Plant(s)

  sp$add_seed()

  expect_that(sp$size, equals(1))
  expect_that(sp$ode_size, equals(p$ode_size))
  expect_that(length(sp$ode_state), equals(p$ode_size))
  expect_that(sp$ode_state, equals(p$ode_state))
  expect_that(sp$ode_rates, equals(rep(NA_real_, p$ode_size)))

  sp$compute_vars_phys(env)
  p$compute_vars_phys(env)
  expect_that(sp$ode_rates, equals(p$ode_rates))
  expect_that(all(sp$ode_rates > 0), is_true())

  pl <- sp$plants
  expect_that(length(pl), equals(1))
  expect_that(class(pl[[1]]), equals(class(p)))
  expect_that(pl[[1]]$ode_state, equals(p$ode_state))
  expect_that(sp$height_max, equals(p$height))
})

test_that("Multiple individuals", {
  h <- 10
  env <- test_environment(h, seed_rain=1.0)
  s <- FFW16_Strategy()
  sp <- FFW16_StochasticSpecies(s)
  n <- 10
  for (i in seq_len(n)) {
    sp$add_seed()
  }

  expect_that(sp$size, equals(n))

  hh <- sort(runif(n, sp$height_max, h), decreasing=TRUE)
  expect_that(sp$heights <- rev(hh), throws_error("must be decreasing"))
  sp$heights <- hh
  expect_that(sp$heights, equals(hh))

  for (i in seq_len(n)) {
    expect_that(sp$plant_at(i)$height, equals(hh[[i]]))
  }

  n_ode <- sp$plant_at(1)$ode_size
  m <- matrix(sp$ode_state, n_ode)
  expect_that(max(m[-1, ]), equals(0))

  nd <- sp$deaths()
  expect_that(nd, equals(0))
  expect_that(sp$size, equals(n))

  m <- matrix(sp$ode_state, n_ode)
  ## Basically set up individual 4 for death:
  i <- 4
  m[2, i] <- 10000
  ## And set another one to something nonzero but trivial:
  j <- 8
  m[2, j] <- .Machine$double.eps
  sp$ode_state <- m
  expect_that(sp$plant_at(i)$mortality_probability, equals(1))
  expect_that(sp$plant_at(j)$mortality_probability, equals(0))
  expect_that(sp$plant_at(j)$mortality, is_more_than(0.0))

  nd <- sp$deaths()
  expect_that(nd, equals(1))
  expect_that(sp$size, equals(n - 1))
  hh2 <- sapply(sp$plants, function(x) x$height)
  expect_that(hh2, equals(hh[-i]))
  expect_that(sp$plant_at(j)$mortality, is_identical_to(0.0))
})


