source("helper-tree.R")

context("Cohort")

s <- new(Strategy)
pars.s <- s$parameters
p <- new(Plant, s)

coh <- new(Cohort, s)

## Probably out-of-order for now, but this checks that initialisation
## happens correctly.
expect_that(coh$ode_values,
            is_identical_to(c(0.0, rep(p$ode_values, 2))))

## Delete to make sure we don't crash on cleanup
rm(coh)
gc()
