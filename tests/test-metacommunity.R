source("helper-tree.R")

context("Metacommunity [Plant]")

p <- new(Parameters)
p$add_strategy(new(Strategy))

sys <- new(Metacommunity, p)
expect_that(sys$size, equals(p$parameters$n_patches))

n <- matrix(0, p$size, sys$size)
expect_that(sys$n_individuals, equals(n))

## Age of system should be zero
expect_that(sys$age, equals(0.0))

## and we should have no ode variables:
expect_that(sys$ode_size, equals(0))

## Put a single individual in the 3rd patch:
idx <- 3
n[idx] <- 1
sys$add_seedlings(n)
expect_that(sys$n_individuals, equals(n))
expect_that(sys$ode_size, equals(sum(n)*3))

expect_that(sys[[idx]]$n_individuals, equals(1)) # 1 individual
expect_that(sys[[idx]]$size, equals(p$size))     # 1 species
expect_that(sys[[idx]][[1]]$size, equals(1))     # 1 plant

y <- new(Plant, p[[1]])$ode_values
expect_that(sys$ode_values, is_identical_to(sys[[idx]]$ode_values))
expect_that(sys$ode_values, is_identical_to(y))

sys$clear()
n[idx] <- 0
expect_that(sys$n_individuals, equals(n))
n[] <- 1
sys$add_seedlings(n)
expect_that(sys$n_individuals, equals(n))

expect_that(sys$ode_values,
            is_identical_to(rep(y, sys$size)))

## Compare derivative calculations for full system and patch.
patch <- new(Patch, p)
patch$add_seedlings(1)
expect_that(sys$ode_values,
            is_identical_to(rep(patch$ode_values, sys$size)))

expect_that(sys$derivs(sys$age, sys$ode_values),
            is_identical_to(rep(patch$derivs(patch$age, patch$ode_values),
                                sys$size)))

expect_that(sys[[idx]]$derivs(sys[[idx]]$age, sys[[idx]]$ode_values),
            is_identical_to(patch$derivs(patch$age, patch$ode_values)))

## TODO: Add a r_set_mass_leaf function that can be used to set leaf
## mass like Species::r_set_mass_leaf and Patch::r_set_mass_leaf.
## Will be easier if we have a get_mass_leaf and if the existing
## functions are modified to take matrices/lists as appropriate.

rm(sys)
gc()
