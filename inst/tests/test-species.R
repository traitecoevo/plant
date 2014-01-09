source("helper-tree.R")
options(error=traceback)

context("Species [Plant]")

s <- new(Strategy)
cmp <- new(Plant, s)

sp <- new(Species, s)
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))
expect_that(sp$height_max, is_identical_to(cmp$height))

sp$add_seeds(1)
expect_that(sp$size, equals(1))
expect_that(sp$n_individuals, equals(1))
expect_that(sp$ode_size, equals(3))

plants <- sp$plants
expect_that(length(plants), equals(1))
expect_that(plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

## Check that access via "[[" works.
expect_that(sp[[1]]$vars_size,
            is_identical_to(cmp$vars_size))
expect_that(sp[[2]]$vars_size,
            throws_error())

## Height to set things to.
h0 <- 10

cmp$height <- h0

expect_that(sp$height <- numeric(0),
            throws_error())
expect_that(sp$height <- c(h0, h0),
            throws_error())
sp$height <- h0

expect_that(sp$plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))
expect_that(sp[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

env <- test.environment(sp$height_max)

sp$compute_vars_phys(env)
cmp$compute_vars_phys(env)

expect_that(sp$plants[[1]]$vars_phys,
            is_identical_to(cmp$vars_phys))

seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(seed$germination_probability(env)))

## Bizarrely, this includes about 100 points; more than I'd have
## thought, and not sure why.  Could be the error in assimilation
## calculation?  Could just be that we're being a bit enthusiastic
## about error checking...
sp$compute_assimilation_fn(env)
tmp <- sp$assimilation_fn
xy <- tmp$xy

## OK, repeat this for Species<CohortDiscrete> and Species<CohortTop>
## and we're done.  For those, remember to tweak the environment times
## and number of indivuduals.
test_that("State get/set works", {
  sp2 <- new(Species, sp$strategy)
  state <- sp2$state
  expect_that(state, is_a("matrix"))
  expect_that(state, equals(matrix(0.0, 3, 0)))

  for (i in 1:4)
    sp2$add_seeds(1)
  sp2$height <- sort(runif(sp2$size), decreasing=TRUE)

  cmp <- rbind(sp2$height, matrix(0, 2, sp2$size))
  expect_that(sp2$state, is_identical_to(cmp))

  tmp <- matrix(sp2$ode_values, nrow=3)
  tmp[2:3,] <- runif(length(tmp[2:3,]))
  sp2$set_ode_values(0, tmp)
  expect_that(sp2$state, is_identical_to(tmp))

  tmp <- tmp + runif(length(tmp))
  tmp[1,] <- sort(tmp[1,], decreasing=TRUE)
  sp2$state <- tmp
  expect_that(sp2$state, is_identical_to(tmp))

  sp2$force_state(tmp[,1:2])
  expect_that(sp2$state, is_identical_to(tmp[,1:2]))

  sp2$force_state(tmp[,2:4])
  expect_that(sp2$state, is_identical_to(tmp[,2:4]))
})

test_that("Size extraction works", {
  sp2 <- new(Species, sp$strategy)
  seed <- sp2$seed

  expect_that(sp2$vars_size, is_a("matrix"))
  expect_that(nrow(sp2$vars_size), equals(length(seed$vars_size)))
  expect_that(ncol(sp2$vars_size), equals(0))
  expect_that(rownames(sp2$vars_size), equals(names(seed$vars_size)))

  ## Add four plants with random sizes.
  for (i in 1:4)
    sp2$add_seeds(1)
  sp2$height <- sort(runif(sp2$size), decreasing=TRUE)

  expect_that(sp2$vars_size, is_a("matrix"))
  expect_that(nrow(sp2$vars_size), equals(length(seed$vars_size)))
  expect_that(ncol(sp2$vars_size), equals(sp2$size))
  expect_that(rownames(sp2$vars_size), equals(names(seed$vars_size)))

  ## Compare the actual values with manually extracting them from the
  ## plants.
  cmp <- sapply(sp2$plants, function(p) p$vars_size)
  expect_that(sp2$vars_size, is_identical_to(cmp))
})

test_that("Physiological variable extraction works", {
  sp2 <- new(Species, sp$strategy)
  seed <- sp2$seed

  expect_that(sp2$vars_phys, is_a("matrix"))
  expect_that(nrow(sp2$vars_phys), equals(length(seed$vars_phys)))
  expect_that(ncol(sp2$vars_phys), equals(0))
  expect_that(rownames(sp2$vars_phys), equals(names(seed$vars_phys)))

  ## Add four plants with random sizes.
  for (i in 1:4)
    sp2$add_seeds(1)
  sp2$height <- sort(runif(sp2$size), decreasing=TRUE)
  sp2$compute_vars_phys(env)

  expect_that(sp2$vars_phys, is_a("matrix"))
  expect_that(nrow(sp2$vars_phys), equals(length(seed$vars_phys)))
  expect_that(ncol(sp2$vars_phys), equals(sp2$size))
  expect_that(rownames(sp2$vars_phys), equals(names(seed$vars_phys)))

  ## Compare the actual values with manually extracting them from the
  ## plants.
  cmp <- sapply(sp2$plants, function(p) p$vars_phys)
  expect_that(sp2$vars_phys, is_identical_to(cmp))
})

## Test approximate plant:

ctrl <- new(Control,
            list(plant_assimilation_approximate_use=TRUE))
s.approx <- new(Strategy)
s.approx$control <- ctrl
sp.approx <- new(Species, s.approx)

sp.approx$add_seeds(1)
sp.approx$height <- sp$height

## Spline is empty on initialisation
expect_that(sp.approx$assimilation_fn$size, equals(0))

sp.approx$compute_vars_phys(env)

## Check that we did compute the spline:
expect_that(sp.approx$assimilation_fn$xy, is_identical_to(tmp$xy))

## The rates should be equal, but not identical, to the fully computed
## rates.
expect_that(sp.approx$ode_rates, equals(sp$ode_rates))
expect_that(identical(sp.approx$ode_rates, sp$ode_rates),
            is_false())

rm(sp)
gc()
