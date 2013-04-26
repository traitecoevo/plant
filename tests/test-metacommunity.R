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
n[3] <- 1
sys$add_plants(n)
expect_that(sys$n_individuals, equals(n))
expect_that(sys$ode_size, equals(sum(n)*3))

rm(sys)
gc()
