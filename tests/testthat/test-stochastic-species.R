context("StochasticSpecies")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("empty", {
  for (x in names(strategy_types)) {
    env <- test_environment(x, 3, seed_rain=1.0)
    e <- environment_types[[x]]
    s <- strategy_types[[x]]()
    sp <- StochasticSpecies(x, e)(s)

    expect_is(sp, sprintf("StochasticSpecies<%s,%s>",x,e))
    expect_equal(sp$size, 0)
    expect_equal(sp$size_plants, 0)

    seed <- sp$seed
    expect_is(seed, "Individual")
    expect_is(seed, sprintf("Individual<%s,%s>",x,e))

    expect_equal(sp$heights, numeric(0))
    expect_equal(sp$height_max, 0.0)
    expect_identical(sp$species, NULL)
    expect_equal(sp$plants, list())
    expect_equal(sp$is_alive, logical())
    expect_equal(sp$compute_competition(0), 0.0)
    expect_equal(sp$ode_size, 0)
    expect_identical(sp$ode_state, numeric(0))
    expect_identical(sp$ode_rates, numeric(0))
  }
})

test_that("Single individual", {
  for (x in names(strategy_types)) {
    env <- test_environment(x, 3, seed_rain=1.0)
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    sp <- StochasticSpecies(x, e)(s)
    p <- Individual(x, e)(s)

    sp$add_seed()

    expect_equal(sp$size, 1)
    expect_equal(sp$size_plants, 1)
    expect_equal(sp$ode_size, p$ode_size)
    expect_equal(length(sp$ode_state), p$ode_size)
    expect_equal(sp$ode_state, p$ode_state)
    expect_equal(sp$ode_rates, rep(NA_real_, p$ode_size))
    expect_equal(sp$is_alive, TRUE)

    sp$compute_rates(env)
    p$compute_rates(env)
    expect_equal(sp$ode_rates, p$ode_rates)

    if (x == "FF16") {
      expect_true(all(sp$ode_rates > 0.0))
    } else if (x == "FF16r") {
      expect_true(all(sp$ode_rates[-3] > 0.0))
      expect_identical(sp$ode_rates[[3]], 0.0)
    }

    pl <- sp$plants
    expect_equal(length(pl), 1)
    expect_equal(class(pl[[1]]), class(p))
    expect_equal(pl[[1]]$ode_state, p$ode_state)
    expect_equal(sp$height_max, p$state("height"))
  }
})

test_that("Multiple individuals", {
  for (x in names(strategy_types)) {
    h <- 10
    env <- test_environment(x, h, seed_rain=1.0)
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    sp <- StochasticSpecies(x, e)(s)
    n <- 10
    for (i in seq_len(n)) {
      sp$add_seed()
    }

    expect_equal(sp$size, n)
    expect_equal(sp$size_plants, n)
    expect_equal(sp$is_alive, rep(TRUE, n))

    hh <- sort(runif(n, sp$height_max, h), decreasing=TRUE)
    expect_error(sp$heights <- rev(hh), "must be decreasing")
    sp$heights <- hh
    expect_equal(sp$heights, hh)

    for (i in seq_len(n)) {
      expect_equal(sp$plant_at(i)$state("height"), hh[[i]])
    }

    n_ode <- sp$plant_at(1)$ode_size
    m <- matrix(sp$ode_state, n_ode)
    expect_identical(max(m[-1, ]), 0.0)

    nd <- sp$deaths()
    expect_equal(nd, 0)
    expect_equal(sp$size, n)
    expect_equal(sp$size_plants, n)

    m <- matrix(sp$ode_state, n_ode)
    ## Basically set up individual 4 for death:
    i <- 4
    m[2, i] <- 10000
    ## And set another one to something nonzero but trivial:
    j <- 8
    m[2, j] <- .Machine$double.eps
    sp$ode_state <- m
    expect_equal(sp$plant_at(i)$mortality_probability, 1)
    expect_equal(sp$plant_at(j)$mortality_probability, 0)
    expect_gt(sp$plant_at(j)$state("mortality"), 0.0)
    m[2, j] <- 0.0 # reset back to original for later comparison

    nd <- sp$deaths()
    expect_equal(nd, 1)
    expect_equal(sp$size, n - 1)
    expect_equal(sp$size_plants, n)
    expect_equal(sp$is_alive, seq_len(n) != i)

    hh2 <- sapply(sp$plants, function(x) x$state("height"))
    ## still the same:
    expect_equal(hh2, hh)
    expect_identical(sp$plant_at(j)$state("mortality"), 0.0)

    m2 <- matrix(sp$ode_state, n_ode)
    expect_equal(ncol(m2), n - 1)
    expect_identical(m2, m[, -i])
  }
})

test_that("establishment probability", {
  for (x in names(strategy_types)) {
    env <- test_environment(x, 3, seed_rain=1.0)
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    sp <- StochasticSpecies(x, e)(s)
    p <- Individual(x, e)(s)

    expect_equal(sp$establishment_probability(env), p$establishment_probability(env))
  }
})
