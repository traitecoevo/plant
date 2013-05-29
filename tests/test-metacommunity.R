source("helper-tree.R")

context("Metacommunity [Plant]")

p <- new(Parameters)
p$add_strategy(new(Strategy))

sys <- new(Metacommunity, p)
expect_that(sys$size, equals(p$parameters$n_patches))

n <- matrix(0, p$size, sys$size)
expect_that(sys$n_individuals, equals(n))

## System time should be zero.
expect_that(sys$time, equals(0.0))

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

expect_that(sys$derivs(sys$time, sys$ode_values),
            is_identical_to(rep(patch$derivs(patch$time, patch$ode_values),
                                sys$size)))

expect_that(sys[[idx]]$derivs(sys[[idx]]$time, sys[[idx]]$ode_values),
            is_identical_to(patch$derivs(patch$time, patch$ode_values)))

expect_that(sys$height,
            equals(rep(list(list(y[[1]])), sys$size)))

## Check dispersal:

## Sanity check -- should fail with bad input
expect_that(sys$disperse(integer(0)),           throws_error())
expect_that(sys$disperse(rep(1, sys$size + 1)), throws_error())

## Compare fit with a Chisq goodness of fit test (may be too harsh?)
n.seeds <- 100
set.seed(1)
m <- replicate(500, sys$disperse(n.seeds)[1,])

## Conservation of seed mass:
expect_that(all(colSums(m) == n.seeds), is_true())

## Goodness of fit.
expected <- n.seeds/sys$size
chisq <- sum((m - expected)^2 / expected)
expect_that(pchisq(chisq, length(m) - 3, lower.tail=FALSE),
            is_greater_than(0.05))


rm(sys)
gc()
