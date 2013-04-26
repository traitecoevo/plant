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
sys$add_plants(n)
expect_that(sys$n_individuals, equals(n))
expect_that(sys$ode_size, equals(sum(n)*3))

patch <- sys[[idx]]
expect_that(patch$n_individuals, equals(1)) # 1 individual
expect_that(patch$size, equals(p$size))     # 1 species
expect_that(patch[[1]]$size, equals(1))     # 1 plant

rm(sys)
gc()
