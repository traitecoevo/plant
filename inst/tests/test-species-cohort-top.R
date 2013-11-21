source("helper-tree.R")

context("Species [CohortTop]")

s <- new(Strategy)
seed <- new(CohortTop, s)
h0 <- seed$height

sp <- new(SpeciesCT, s)

## Note that this initialises with a special base cohort always.
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))
expect_that(sp$height_max, is_identical_to(seed$height))

env <- test.environment(seed$height * 10, seed.rain=1.0)

## Causes initial conditions to be estimated:
sp$compute_vars_phys(env)
seed$compute_initial_conditions(env)

sp$add_seeds(1)

expect_that(sp[[1]]$vars_phys,
            is_identical_to(seed$vars_phys))
expect_that(sp$ode_values,
            is_identical_to(seed$ode_values))

## Internal and test seed report same values:
expect_that(sp$seed$vars_phys,
            is_identical_to(seed$vars_phys))
expect_that(sp$seed$ode_values,
            is_identical_to(seed$ode_values))
## Finished with that seed now, so delete.
rm(seed)

plant.seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(plant.seed$germination_probability(env)))

expect_that(sp$leaf_area_above(0), equals(0))

sp$height <- 1
h <- 0
x <- c(sp$seed$height, sp$height)
y <- c(sp$seed$leaf_area_above(h), sp[[1]]$leaf_area_above(h))

expect_that(sp$leaf_area_above(h),
            is_identical_to(trapezium(x, y)))

## Better tests: I want cases where:
## 1. empty: throws error

sp$clear()

## Re-set up the initial conditions
sp$compute_vars_phys(env)

expect_that(sp[[1]], throws_error())

## TODO: What to do about negative height?

## 1: empty species (no cohorts) has no leaf area above any height:
test_that("Empty species has no leaf area", {
  expect_that(sp$leaf_area_above(0), equals(0))
  expect_that(sp$leaf_area_above(10), equals(0))
  expect_that(sp$leaf_area_above(Inf), equals(0))
})

## 2: Cohort up against boundary has no leaf area:
sp$add_seeds(1)
test_that("Species with only boundary cohort no leaf area", {
  expect_that(sp$leaf_area_above(0), equals(0))
  expect_that(sp$leaf_area_above(10), equals(0))
  expect_that(sp$leaf_area_above(Inf), equals(0))
})

## 3: Single cohort; one round of trapezium:
h.top <- h0 * 4
sp$height <- h.top
cmp <- function(h, sp) {
  x <- c(sp$height, sp$seed$height)
  y <- c(sapply(sp$plants, function(p) p$leaf_area_above(h)),
         sp$seed$leaf_area_above(h))
  trapezium(rev(x), rev(y))
}

test_that("Leaf area sensible with one cohort", {
  ## At base and top
  expect_that(sp$leaf_area_above(0), equals(cmp(0, sp)))
  expect_that(sp$leaf_area_above(h.top), equals(0))
  ## Part way up (and above bottom seed boundary condition)
  expect_that(sp$leaf_area_above(h.top * .5),
              equals(cmp(h.top * .5, sp)))
})

sp$add_seeds(1)
sp$height <- h.top * c(1, .6)

test_that("Leaf area sensible with two cohorts", {
  ## At base and top
  expect_that(sp$leaf_area_above(0), equals(cmp(0, sp)))
  expect_that(sp$leaf_area_above(h.top), equals(0))
  ## Part way up (below bottom cohort, above boundarty condition)
  expect_that(sp$leaf_area_above(h.top * .5),
              equals(cmp(h.top * .5, sp)))
  ## Within the top pair (excluding the seed)
  expect_that(sp$leaf_area_above(h.top * .8),
              equals(cmp(h.top * .8, sp)))
})

sp$add_seeds(1)
sp$height <- h.top * c(1, .75, .6)
test_that("Leaf area sensible with two cohorts", {
  ## At base and top
  expect_that(sp$leaf_area_above(0), equals(cmp(0, sp)))
  expect_that(sp$leaf_area_above(h.top), equals(0))
  ## Part way up (below bottom cohort, above boundarty condition)
  expect_that(sp$leaf_area_above(h.top * .5),
              equals(cmp(h.top * .5, sp)))
  ## Within the top pair (excluding the seed)
  expect_that(sp$leaf_area_above(h.top * .8),
              equals(cmp(h.top * .8, sp)))
})

test_that("Total leaf area calculations are correct", {
  cmp.leaf.area <- sapply(seq_len(sp$size), function(i) sp[[i]]$leaf_area)
  expect_that(sp$leaf_area,
              is_identical_to(cmp.leaf.area))
})

hh <- sp$height
hh[1] <- 10
sp$height <- hh
env <- test.environment(max(hh), seed.rain=1.0)
sp$compute_assimilation_spline(env)
tmp <- sp$assimilation_spline

h.sp <- tmp$x
h.sp.mid <- (h.sp[-1] + h.sp[-length(h.sp)]) / 2
a.sp <- tmp$y

pl <- new(Plant, s)
a.pl <- sapply(h.sp, function(h) pl$assimilation_given_height(h, env))
expect_that(a.pl, is_identical_to(a.sp))

a.sp.mid <- tmp$eval(h.sp.mid)
a.pl.mid <- sapply(h.sp.mid, function(h)
                   pl$assimilation_given_height(h, env))

## OK; this is surprising because we're seeing step changes even
## though we apparently never subdvide.  The oscillation at the end is
## probably OK given the general scaling of the solution though.
## These might just be where the level of detail has changed though.
## It does look like the error estimate might be wrong though -- we're
## going bonkers at the bottom end; I should be happy if *either*
## absolute or relative error are satisfied, but it looks like I might
## be driving *both* down?
## plot(a.sp.mid, a.pl.mid)
## plot(a.sp.mid, a.pl.mid - a.sp.mid)
## plot(h.sp.mid, a.pl.mid - a.sp.mid)

test_that("Approximately computed assimilation matches full computation", {
  ctrl.full <- new(Control,
                   list(cohort_gradient_richardson=TRUE))
  s.full <- new(Strategy)
  s.full$control <- ctrl.full
  sp.full <- new(SpeciesCT, s.full)

  ctrl.approx <- new(Control,
                     list(plant_assimilation_approximate_use=TRUE))
  s.approx <- new(Strategy)
  s.approx$control <- ctrl.approx
  sp.approx <- new(SpeciesCT, s.approx)

  sp.full$compute_vars_phys(env)
  sp.approx$compute_vars_phys(env)

  sp.approx$add_seeds(3)
  sp.full$add_seeds(3)

  sp.approx$height <- sp$height
  sp.full$height <- sp$height

  sp.approx$compute_vars_phys(env)
  sp.full$compute_vars_phys(env)

  expect_that(sp.approx$ode_values, equals(sp.full$ode_values))

  ## The rates should be equal, but not identical, to the fully computed
  ## rates.
  ## OK, these are actually quite different, but that's because the ODE
  ## values are also quite different, including for the seed.
  r.approx <- matrix(sp.approx$ode_rates, 4)
  r.full   <- matrix(sp.full$ode_rates,   4)

  expect_that(r.approx[1:3,], equals(r.full[1:3,], tolerance=2e-7))
  expect_that(r.approx[4,],   equals(r.full[4,],   tolerance=0.0008))
  expect_that(identical(r.full, r.approx), is_false())
})

test_that("State get/set works", {
  set.seed(1)
  sp2 <- new(SpeciesCT, sp$strategy)
  env$time <- 0
  sp2$compute_vars_phys(env)

  seed <- new(CohortTop, s)
  seed$compute_initial_conditions(env)

  state <- sp2$state
  expect_that(state, is_a("matrix"))
  expect_that(state, is_identical_to(cbind(seed$state)))

  ## First few just check that the expected behaviour with the initial
  ## conditions work.
  sp2$add_seeds(1)
  sp2$height <- h.top * 0.5
  state <- cbind(state, state)
  state[1,1] <- sp2$height
  expect_that(sp2$state, is_identical_to(state))

  env$time <- 10
  sp2$compute_vars_phys(env)
  seed$compute_initial_conditions(env)
  state[,2] <- seed$state
  expect_that(sp2$state, is_identical_to(state))

  sp2$add_seeds(1)
  env$time <- 20
  sp2$compute_vars_phys(env)
  seed$compute_initial_conditions(env)
  state <- cbind(state, seed$state)
  expect_that(sp2$state, is_identical_to(state))

  ## Then set some values and check that they stick:
  state[1,1:2] <- sort(runif(2), decreasing=TRUE)
  state[2:5,] <- runif(4 * ncol(state))
  sp2$state <- state
  expect_that(sp2$state, is_identical_to(state))

  ## Check that we can't set invalid states:
  expect_that(sp2$state <- state[,1:2],             throws_error())
  expect_that(sp2$state <- cbind(state[,1], state), throws_error())
  expect_that(sp2$state <- state[1:4,],             throws_error())
  expect_that(sp2$state <- rbind(state[1,], state), throws_error())

  ## Then look at forcing state:
  state <- cbind(state[,1], state)
  state[1,-ncol(state)] <- sort(runif(ncol(state)-1), decreasing=TRUE)
  sp2$force_state(state)
  expect_that(sp2$state,      is_identical_to(state))
  expect_that(sp2$seed$state, is_identical_to(state[,ncol(state)]))
})

rm(sp)
gc()
