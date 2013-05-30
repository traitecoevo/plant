source("helper-tree.R")

context("Patch [CohortTop]")

p <- new(Parameters)
p$add_strategy(new(Strategy))

## An individual CohortTop that will be the same in terms of strategy
## (and initial mass).
cmp <- new(CohortTop, p[[1]])

patch.p <- new(Patch,  p)
patch.c <- new(PatchCohortTop, p)

expect_that(patch.p$height_max, is_identical_to(cmp$height))
expect_that(patch.c$height_max, is_identical_to(cmp$height))

expect_that(patch.p$size, equals(1))
expect_that(patch.p$ode_size, equals(0))

r <- pi/2
patch.c$set_seed_rain(seed_rain(r))
expect_that(patch.c$environment$seed_rain$seed_rain,
            is_identical_to(r))

## This should not work, because we can't set environment from within
## patch.
expect_that(patch.c$environment$seed_rain <- seed_rain(r * 2),
            throws_error())

## This is required right now, I *think*.
## patch.c$compute_vars_phys()

patch.c$add_seedling(1)
patch.p$add_seedling(1)
expect_that(patch.c$ode_size, equals(4))

## Now that we've got started, we should not be able to set the seed
## rain:
expect_that(patch.c$set_seed_rain(seed_rain(1.0)),
            throws_error())

cmp$compute_initial_conditions(patch.c$environment)

expect_that(patch.c$ode_values,
            is_identical_to(cmp$ode_values))
expect_that(patch.c$ode_rates,
            is_identical_to(cmp$ode_rates))

## This gives a "requested rain out of bounds", which is probably the
## cause of the EBT error -- nice.
##
## So, update to get the seed rain initialised, and we're probably
## sorted.
##
## I'm still a bit confused though, as to why this crashes the ebt...
## patch.c$set_ode_values(0, y)
