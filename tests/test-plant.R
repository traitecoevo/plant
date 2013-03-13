source("helper-tree.R")

context("Plant")

cmp <- make.falster()
pars.cmp <- cmp$get_parameters()

s <- new(Strategy)
core <- list(lma=0.1978791, hmat=16.5958691, rho=608, s=3.8e-5)
s$set_params(core)
pars.s <- s$get_params()

## Expect that all parameters in the R version are found in the C++
## version.
expect_that(all(names(pars.cmp) %in% names(pars.s)), is_true())

## And v.v.
expect_that(all(names(pars.s) %in% names(pars.cmp)), is_true())

## The C++ version should have no NA values by this point.
expect_that(any(sapply(pars.s, is.na)), is_false())

## But the R version may have a few:
pars.na <- sort(names(pars.cmp)[sapply(pars.cmp, is.na)])
expect_that(pars.na,
            is_identical_to(c("c_d0", "c_d1", "c_d2", "c_d3", 
                              "c_s0", "Pi_0", "s")))

## And demand that all parameters agree.
pars.ok <- setdiff(names(pars.cmp), pars.na)
expect_that(pars.s[pars.ok], is_identical_to(pars.cmp[pars.ok]))

## Now that that's OK, make a plant with our strategy
p <- new(Plant, s)
## ...so the parameters should agree exactly.
expect_that(p$parameters, is_identical_to(pars.s))

## Set the leaf mass to something (here pi)
p$set_mass_leaf(pi)

size.p <- p$vars_size

## Daniel's is done in terms of height, so start there:
h <- size.p[["height"]]
a <- size.p[["leaf_area"]]

expect_that(cmp$LeafMass(cmp$traits$lma, cmp$LeafArea(h)),
            equals(pi))

expect_that(size.p[["leaf_area"]],
            equals(cmp$LeafArea(h)))
expect_that(size.p[["height"]],
            equals(cmp$Height(a)))
expect_that(size.p[["mass_leaf"]], equals(pi))
expect_that(size.p[["mass_sapwood"]],
            equals(cmp$SapwoodMass(cmp$traits$rho, a, h)))
## This one differs slightly.
expect_that(size.p[["mass_heartwood"]],
            equals(cmp$HeartwoodMass(cmp$traits$rho, a)))
expect_that(size.p[["mass_root"]],
            equals(cmp$RootMass(a)))
expect_that(size.p[["mass_total"]],
            equals(cmp$TotalMass(cmp$traits, a)))

## Delete the plant -- should not crash.
rm(p)
gc()

