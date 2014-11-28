source("helper-tree.R")
options(error=traceback)

context("Species [CohortDiscrete]")

s <- new(Strategy)
cmp <- new(Plant, s)
ode_size_plant <- cmp$ode_size

sp <- new(SpeciesC, s)
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))
expect_that(sp$height_max, is_identical_to(cmp$height))

sp$add_seeds(2)
expect_that(sp$size, equals(1))
expect_that(sp$n_individuals, equals(2))
expect_that(sp$ode_size, equals(ode_size_plant))

expect_that(sp$size, equals(1))
expect_that(sp[[1]]$vars_size,
            is_identical_to(cmp$vars_size))
expect_that(length(sp$plants), equals(1))

h0 <- 10
cmp$height <- h0

expect_that(sp$height <- numeric(0),
            throws_error())
expect_that(sp$height <- c(h0, h0),
            throws_error())
sp$height <- h0

expect_that(sp[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

env <- test.environment(sp$height_max)

sp$compute_vars_phys(env)
cmp$compute_vars_phys(env)

expect_that(sp[[1]]$vars_phys,
            is_identical_to(cmp$vars_phys))

seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(seed$germination_probability(env)))

test_that("State get/set works", {
  state_size <- ode_size_plant+1  # Not sure why - Rich?
  set.seed(1)
  sp2 <- new(SpeciesC, sp$strategy)
  state <- sp2$state
  expect_that(state, is_a("matrix"))
  expect_that(state, equals(matrix(0.0, state_size, 0)))

  sp2$add_seeds(1)
  state <- sp2$state
  expect_that(state,
              is_identical_to(matrix(c(seed$height, 0, 0, 0, 0, 1), state_size, 1)))

  sp2$add_seeds(2)
  expect_that(sp2$state,
              is_identical_to(cbind(state, c(seed$height, 0, 0, 0, 0, 2))))
  state <- sp2$state

  sp2$height <- state[1,] <- sort(runif(2), decreasing=TRUE)
  expect_that(sp2$state, is_identical_to(state))

  state[2:3,] <- runif(4)
  sp2$set_ode_values(0, state[1:ode_size_plant,])
  expect_that(sp2$state, is_identical_to(state))

  ## Now, set to a brand new state:
  state[1,] <- sort(runif(2), decreasing=TRUE) * (2*ode_size_plant)
  state[2:3,] <- runif(4)
  state[4,] <- c(3, 6)
  sp2$state <- state
  expect_that(sp2$state, is_identical_to(state))

  state <- cbind(state, c(seed$height, 0, 0, 0, 0, 2))
  expect_that(sp2$state <- state,                throws_error())
  expect_that(sp2$state <- state[,1,drop=FALSE], throws_error())
  expect_that(sp2$state <- state[1:3,],          throws_error())
  expect_that(sp2$state <- state[c(1:4,1),],     throws_error())

  sp2$force_state(state)
  expect_that(sp2$state, is_identical_to(state))
})

rm(sp)
gc()
