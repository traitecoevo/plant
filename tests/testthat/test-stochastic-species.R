
## TODO: Remove ["FFW16"] to test this with all types.
## But this triggers error on `all(sp$ode_rates > 0)`
strategy_types <- get_list_of_strategy_types()["FFW16"]

for (x in names(strategy_types)) {

  context(sprintf("StochasticSpecies-%s",x))

  test_that("empty", {
    env <- test_environment(3, seed_rain=1.0)
    s <- strategy_types[[x]]()
    sp <- StochasticSpecies(x)(s)

    expect_that(sp, is_a(sprintf("StochasticSpecies<%s>",x)))
    expect_that(sp$size, equals(0))
    expect_that(sp$size_plants, equals(0))

    seed <- sp$seed
    expect_that(seed, is_a("Plant"))
    expect_that(seed, is_a(sprintf("Plant<%s>",x)))

    expect_that(sp$heights, equals(numeric(0)))
    expect_that(sp$height_max, equals(0.0))
    expect_that(sp$species, is_identical_to(NULL))
    expect_that(sp$plants, equals(list()))
    expect_that(sp$is_alive, equals(logical()))
    expect_that(sp$area_leaf_above(0), equals(0.0))
    expect_that(sp$ode_size, equals(0))
    expect_that(sp$ode_state, is_identical_to(numeric(0)))
    expect_that(sp$ode_rates, is_identical_to(numeric(0)))
  })

  test_that("Single individual", {
    env <- test_environment(3, seed_rain=1.0)
    s <- strategy_types[[x]]()
    sp <- StochasticSpecies(x)(s)
    p <- Plant(x)(s)

    sp$add_seed()

    expect_that(sp$size, equals(1))
    expect_that(sp$size_plants, equals(1))
    expect_that(sp$ode_size, equals(p$ode_size))
    expect_that(length(sp$ode_state), equals(p$ode_size))
    expect_that(sp$ode_state, equals(p$ode_state))
    expect_that(sp$ode_rates, equals(rep(NA_real_, p$ode_size)))
    expect_that(sp$is_alive, equals(TRUE))

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
    s <- strategy_types[[x]]()
    sp <- StochasticSpecies(x)(s)
    n <- 10
    for (i in seq_len(n)) {
      sp$add_seed()
    }

    expect_that(sp$size, equals(n))
    expect_that(sp$size_plants, equals(n))
    expect_that(sp$is_alive, equals(rep(TRUE, n)))

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
    expect_that(sp$size_plants, equals(n))

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
    m[2, j] <- 0.0 # reset back to original for later comparison

    nd <- sp$deaths()
    expect_that(nd, equals(1))
    expect_that(sp$size, equals(n - 1))
    expect_that(sp$size_plants, equals(n))
    expect_that(sp$is_alive, equals(seq_len(n) != i))

    hh2 <- sapply(sp$plants, function(x) x$height)
    ## still the same:
    expect_that(hh2, equals(hh))
    expect_that(sp$plant_at(j)$mortality, is_identical_to(0.0))

    m2 <- matrix(sp$ode_state, n_ode)
    expect_that(ncol(m2), equals(n - 1))
    expect_that(m2, is_identical_to(m[, -i]))
  })

  test_that("germination probability", {
    env <- test_environment(3, seed_rain=1.0)
    s <- strategy_types[[x]]()
    sp <- StochasticSpecies(x)(s)
    p <- Plant(x)(s)

    expect_that(sp$germination_probability(env),
                equals(p$germination_probability(env)))
  })
}