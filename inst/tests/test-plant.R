source("helper-tree.R")

context("Plant")

cmp <- make.reference()
pars.cmp <- cmp$get_parameters()

s <- new(Strategy)
pars.s <- s$parameters

## Expect that all parameters in the R version are found in the C++
## version.
expect_that(all(names(pars.cmp) %in% names(pars.s)),
            is_true())

## And v.v.
expect_that(all(names(pars.s) %in% names(pars.cmp)),
            is_true())

## The C++ version should have no NA values by this point.
expect_that(any(sapply(pars.s, is.na)),
            is_false())

## And neither should the R version.
expect_that(any(sapply(pars.cmp, is.na)),
            is_false())

## And demand that all parameters agree.
expect_that(pars.s[names(pars.cmp)],
            is_identical_to(pars.cmp))

## TODO: This should fail, but currently does not.
## Not sure which constructor is being called, or how.
## expect_that(new(Plant),
##             throws_error())
## It's not calling either of the defined constructors, so must be
## using a defalt-defined no-args constructor.  Surprised though.
## Could this be a bug in Rcpp?

## Now that that's OK, make a plant with our strategy
p <- new(Plant, s)
## ...so the parameters should agree exactly.
expect_that(p$strategy$parameters,
            is_identical_to(pars.s))

## Set the height to something (here 10)
h0 <- 10
p$height <- h0

size.p <- p$vars_size

expect_that(size.p[["height"]],
            equals(h0))
expect_that(size.p[["leaf_area"]],
            equals(cmp$LeafArea(h0)))
expect_that(size.p[["mass_leaf"]],
            equals(cmp$LeafMass(cmp$traits$lma, cmp$LeafArea(h0))))
expect_that(size.p[["mass_sapwood"]],
            equals(cmp$SapwoodMass(cmp$traits$rho, cmp$LeafArea(h0), h0)))
expect_that(size.p[["mass_bark"]],
            equals(cmp$BarkMass(cmp$traits$rho, cmp$LeafArea(h0), h0)))
expect_that(size.p[["mass_root"]],
            equals(cmp$RootMass(cmp$LeafArea(h0))))
expect_that(size.p[["mass_live"]],
            equals(cmp$LiveMass(cmp$traits, cmp$LeafArea(h0))))
expect_that(size.p[["area_bark"]],
            equals(cmp$bark_area(h0)))
expect_that(size.p[["area_sapwood"]],
            equals(cmp$sapwood_area(h0)))
expect_that(size.p[["area_basal"]],
            equals(cmp$basal_area(h0)))
expect_that(p$height, is_identical_to(size.p[["height"]]))
expect_that(p$leaf_area, is_identical_to(size.p[["leaf_area"]]))

## check heartwood area function

## Expect zero unless it has been set otherwise
expect_that(size.p[["area_heartwood"]],
            equals(0))
HA0 <- 1E-3
p$heartwood_area <- HA0
size.p <- p$vars_size
expect_that(size.p[["area_heartwood"]],
            equals(HA0))
expect_that(size.p[["area_basal"]],
            equals(cmp$basal_area(h0) + HA0))
# set heartwood back at zero for subsequent tests
p$heartwood_area <- 0

## Light environment

env <- test.environment(h0)
light.env <- attr(env, "light.env") # underlying function

## The R model computes A_lf * leaf_area * Y * c_bio, wheras we just
## compute A_lf; will have to correct some numbers.
cmp.const <- pars.s$Y * pars.s$c_bio

## Compute the physiological variables and retreive them.
p$compute_vars_phys(env)
p.phys <- p$vars_phys
p.growth_decomp <- p$vars_growth_decomp

## 1. Assimilation:
cmp.assimilation.plant <- cmp$assimilation.plant(h0, light.env)
expect_that(p.phys[["assimilation"]],
            equals(cmp.assimilation.plant / cmp.const))

## 2. Respiration:
cmp.respiration <- cmp$respiration.given.height(cmp$traits, h0)
expect_that(p.phys[["respiration"]],
            equals(cmp.respiration / cmp.const))

## 3. Turnover:
cmp.turnover <- cmp$turnover.given.height(cmp$traits, h0)
expect_that(p.phys[["turnover"]],
            equals(cmp.turnover))

## 4. Net production:
cmp.net.production <- cmp$net.production(cmp$traits, h0, light.env)
expect_that(p.phys[["net_production"]],
            equals(cmp.net.production, tolerance=1e-7))

## 5. Reproduction fraction
cmp.reproduction.fraction <-
  cmp$ReproductiveAllocation(cmp$traits$hmat,h0)
expect_that(p.phys[["reproduction_fraction"]],
            equals(cmp.reproduction.fraction))

## 6. Fecundity rate
cmp.fecundity.rate <- cmp$fecundity.rate(cmp$traits, h0, light.env)
expect_that(p.phys[["fecundity_rate"]],
            equals(cmp.fecundity.rate, tolerance=1e-7))

## 7. Fraction of whole plant (mass) growth that is leaf.
cmp.leaf.fraction <- cmp$leaf.fraction(cmp$traits, h0)
expect_that(p.phys[["leaf_fraction"]],
            equals(cmp.leaf.fraction))

## 8. Growth rate for height
cmp.height.growth.rate <- cmp$height.growth.rate(cmp$traits, h0, light.env)
expect_that(p.phys[["height_growth_rate"]],
            equals(cmp.height.growth.rate, tolerance=1e-7))

cmp.height.growth.rate <-
  cmp$height.growth.rate.via.mass.leaf(cmp$traits, h0, light.env)
expect_that(p.phys[["height_growth_rate"]],
            equals(cmp.height.growth.rate, tolerance=1e-7))

## 9. Mortality rate
cmp.mortality.rate <- cmp$mortality.rate(cmp$traits, h0, light.env)
expect_that(p.phys[["mortality_rate"]],
            equals(cmp.mortality.rate))

## 10. Archietcural layout
cmp.dheight_dleaf_area <- cmp$dHdA(cmp$LeafArea(h0))
expect_that(p.growth_decomp[["dheight_dleaf_area"]],
            equals(cmp.dheight_dleaf_area))

## 11. Sapwood mass per leaf mass
cmp.dmass_sapwood_dmass_leaf <- cmp$sapwood.per.leaf.mass(cmp$traits, h0)
expect_that(p.growth_decomp[["dmass_sapwood_dmass_leaf"]],
            equals(cmp.dmass_sapwood_dmass_leaf))

## 12. Bark mass per leaf mass
cmp.dmass_bark_dmass_leaf <- cmp$bark.per.leaf.mass(cmp$traits, h0)
expect_that(p.growth_decomp[["dmass_bark_dmass_leaf"]],
            equals(cmp.dmass_bark_dmass_leaf))

## 12. Root mass per leaf mass
cmp.dmass_root_dmass_leaf <- cmp$root.per.leaf.mass(cmp$traits, h0)
expect_that(p.growth_decomp[["dmass_root_dmass_leaf"]],
            equals(cmp.dmass_root_dmass_leaf))

## 13. Leaf area growth rate
cmp.dleaf_area_dt <- cmp$dleaf_area_dt(cmp$traits, h0, light.env)
expect_that(p.growth_decomp[["dleaf_area_dt"]],
            equals(cmp.dleaf_area_dt, tolerance=1e-7))

## 14. sapwood area growth rate
cmp.dsapwood_area_dt <- cmp$dsapwood_area_dt(cmp$traits, h0, light.env)
expect_that(p.growth_decomp[["dsapwood_area_dt"]],
            equals(cmp.dsapwood_area_dt, tolerance=1e-7))

## 15. bark area growth rate
cmp.dbark_area_dt <- cmp$dbark_area_dt(cmp$traits, h0, light.env)
expect_that(p.growth_decomp[["dbark_area_dt"]],
            equals(cmp.dbark_area_dt, tolerance=1e-7))

## 16. heartwood area growth rate
cmp.dheartwood_area_dt <- cmp$dheartwood_area_dt(cmp$traits, h0, light.env)
expect_that(p.growth_decomp[["dheartwood_area_dt"]],
            equals(cmp.dheartwood_area_dt, tolerance=1e-7))

## 17. basal area growth rate
cmp.dbasal_area_dt <- cmp$dbasal_area_dt(cmp$traits, h0, light.env)
expect_that(p.growth_decomp[["dbasal_area_dt"]],
            equals(cmp.dbasal_area_dt, tolerance=1e-7))

## 18. change in basal diam per basal area
cmp.dbasal_diam_dbasal_area <- cmp$dbasal_diam_dbasal_area(cmp$basal_area(h0))
expect_that(p.growth_decomp[["dbasal_diam_dbasal_area"]],
            equals(cmp.dbasal_diam_dbasal_area, tolerance=1e-7))

## 18. basal diam growth rate
cmp.dbasal_diam_dt <- cmp$dbasal_diam_dt(cmp$traits, h0, light.env)
expect_that(p.growth_decomp[["dbasal_diam_dt"]],
            equals(cmp.dbasal_diam_dt, tolerance=1e-7))

## Check that height decomposition multiplies out to give right answer
expect_that(p.growth_decomp[["height_growth_rate"]],
            equals(
              prod(p.growth_decomp[c("dheight_dleaf_area","dleaf_area_dleaf_mass","leaf_fraction","growth_fraction","net_production")]), tolerance=1e-7))

## Seed stuff:
seed <- new(Plant, s)

## Check that our root-finding succeeded and the leaf mass is correct:
expect_that(seed$vars_size[["mass_live"]],
            equals(pars.s$s, tolerance=1e-7))

## Check that the height at birth is correct.  These answers are
## actually quite different, which could come from the root finding?
expect_that(seed$height,
            equals(cmp$height.at.birth(cmp$traits), tolerance=1e-4))

## Then, check the germination probabilities in the current light
## environment:
expect_that(seed$germination_probability(env),
            equals(cmp$germination.probability(cmp$traits, light.env),
                   tolerance=1e-5))

## Check ode state get/set:
test_that("ODE state has known order", {
  vals <- p$ode_values
  expect_that(length(vals), equals(5))
  expect_that(vals[[1]], is_identical_to(p$height))

  p2 <- new(Plant, p$strategy)
  vals[2:3] <- runif(2)
  p2$set_ode_values(NA, vals)
  expect_that(p2$ode_values, is_identical_to(vals))
  expect_that(vals[[1]], is_identical_to(p2$height))
  expect_that(vals[[2]], equals(-log(p2$survival_probability)))
})

test_that("System state get/set works", {
  expect_that(p$state_size, equals(5))
  vals <- p$state
  expect_that(vals, is_identical_to(p$ode_values))

  vals2 <- vals + runif(5)
  p2 <- new(Plant, p$strategy)
  p2$state <- vals2
  expect_that(p2$state, is_identical_to(vals2))
  expect_that(p2$state, is_identical_to(p2$ode_values))

  expect_that(p2$state <- c(vals, 1), throws_error())
  expect_that(p2$state <- vals[1:2],  throws_error())
})

## Check with different control parameters.
## TODO: This is really awkward.
ctrl <- s$control
ctrl$set_parameters(list(plant_assimilation_over_distribution=TRUE))
s$control <- ctrl
p2 <- new(Plant, s)
p2$height <- p$height

p2$compute_vars_phys(env)
p2.phys <- p2$vars_phys
test_that("Assimilation over distribution works", {
  expect_that(p2.phys, equals(p.phys))
  ## Will be simplifable by devtools::not()
  expect_that(identical(p2.phys, p.phys), is_false())
})

## Grow the plant in a constant environment
derivs <- function(t, y, pars) {
  plant <- pars$plant
  light.env <- pars$light.env

  plant$set_ode_values(t, y)
  plant$compute_vars_phys(light.env)
  plant$ode_rates
}

## Make a bigger light environment, so the plant has room to grow
env2 <- test.environment(pars.s$hmat * 1.2, 300)
p$compute_vars_phys(env2)
p.phys <- p$vars_phys
p.growth_decomp <- p$vars_growth_decomp

## Check the derivative calculations are correct
t <- 0.0 # arbitrary, ignored
y <- c(h0, 0, 0, 0,0)
pars.derivs <- list(plant=p, light.env=env2)
tmp <- derivs(t, y, pars.derivs)
p.derivs <- c(p.phys[c("height_growth_rate",
                     "mortality_rate", "fecundity_rate")],  p.growth_decomp[c("dheartwood_area_dt","dheartwood_mass_dt")] )
expect_that(tmp,
            equals(unname(p.derivs), tolerance=2e-8))

## and again
t0 <- 0.0
obj <- new(OdeR, derivs, new.env(), pars.derivs)
expect_that(obj$derivs(t0, y),
            equals(unname(p.derivs), tolerance=2e-8))

## Then run for 15 years
tt <- seq(0, 15, length=51)

library(deSolve)
derivs.d <- function(...)
  list(derivs(...))
cmp.run <- unname(rk(y, tt, derivs.d, pars.derivs,
                     method=rkMethod("rk45ck"), hini=1/1000, rtol=1,
                     atol=1)[-1,-1])

expect_that(t(obj$run(tt, y)),
            equals(cmp.run, tolerance=1e-7))

## Then, test the one-shot assimilation function.
test_that("One shot assimilation function works", {
  s$control <- new(Control)
  p.cmp <- new(Plant, s)

  set.seed(1)
  hh <- runif(10, 1, 10)
  aa <- sapply(hh, function(h) p.cmp$assimilation_given_height(h, env))
  aa.R <- sapply(hh, function(h) cmp$assimilation.plant(h, light.env))/
    cmp.const

  foo <- function(h, p, env) {
    p$height <- h
    p$compute_vars_phys(env)
    p$vars_phys[["assimilation"]]
  }
  expect_that(aa, is_identical_to(sapply(hh, foo, p, env)))
  expect_that(aa, equals(aa.R))
})

test_that("Approximate assimilation works", {
  ctrl <- new(Control)
  ctrl$set_parameters(list(plant_assimilation_approximate_use=TRUE))
  s2 <- new(Strategy)
  s2$control <- ctrl

  assim <- compute_assimilation_fn(s2, seed$height, h0, env)
  expect_that(assim$max, is_identical_to(h0))

  ## Then set this within the strategy so that it would be available to
  ## a plant.
  s2$assimilation_fn <- assim
  p2 <- new(Plant, s2)

  p2$height <- p$height
  p2$heartwood_area <-  p$heartwood_area
  p2$heartwood_mass <-  p$heartwood_mass

  p$compute_vars_phys(env)
  p2$compute_vars_phys(env)

  expect_that(p2$vars_size, is_identical_to(p$vars_size))

  ## Expect that the physiological variables will be similar, but not
  ## exactly the same.
  expect_that(p2$vars_phys, equals(p$vars_phys))
  expect_that(identical(p2$vars_phys, equals(p$vars_phys)), is_false())
  expect_that(p2$ode_rates, equals(p$ode_rates))
})

test_that("Non-adaptive assimilation integration works", {
  # A new Strategy object that is identical to previous one, but uses
  # non-adaptive integration:
  ctrl <- new(Control, list(plant_assimilation_adaptive=FALSE))
  s.f <- s$copy()
  s.f$control <- ctrl

  p.a <- new(Plant, s)
  p.f <- new(Plant, s.f)

  ## Make the plants more interesting size:
  p.a$height <- 10
  p.f$height <- 10

  p.a$compute_vars_phys(env)
  p.f$compute_vars_phys(env)

  expect_that(p.a$strategy$integrator$is_adaptive, is_true())
  expect_that(p.f$strategy$integrator$is_adaptive, is_false())
})

## Delete the plant -- should not crash.
rm(p)
gc()

## Test copy and assignment on standalone and non-standalone plants
expect_that(test_plant(s, FALSE, FALSE), is_true())
expect_that(test_plant(s, TRUE,  FALSE), is_true())
expect_that(test_plant(s, FALSE, TRUE),  is_true())
expect_that(test_plant(s, TRUE,  TRUE),  is_true())
