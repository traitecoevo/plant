source("helper-tree.R")

context("Plant")

cmp <- make.falster()
pars.cmp <- cmp$get_parameters()

s <- new(Strategy)
pars.s <- s$get_parameters()

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
expect_that(p$parameters,
            is_identical_to(pars.s))

## Expect that any two plants differ in their "name".
expect_that(p$name != new(Plant, s)$name,
            is_true())

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
expect_that(size.p[["mass_bark"]],
            equals(cmp$BarkMass(cmp$traits$rho, a, h)))
expect_that(size.p[["mass_heartwood"]],
            equals(cmp$HeartwoodMass(cmp$traits$rho, a)))
expect_that(size.p[["mass_root"]],
            equals(cmp$RootMass(a)))
expect_that(size.p[["mass_total"]],
            equals(cmp$TotalMass(cmp$traits, a)))

## Instantaneous photosynthesis:
xx <- seq(0, 1, length=101)
yy <- sapply(xx, p$assimilation_leaf)
zz <- cmp$Assim(a, xx)

## The R model computes A_lf * leaf_area * Y * c_bio, wheras we just
## compute A_lf
cmp.const <- pars.s$Y * pars.s$c_bio
expect_that(yy * cmp.const * a,
            equals(zz))

## Make a pretend light environment over the plant height, slightly
## concave up, whatever.
hh <- seq(0, h, length=101)
light.env <- function(x)
  exp(x/(h*2)) - 1 + (1 - (exp(.5) - 1))/2
ee <- light.env(hh)

cmp.assimilation.leaf <- cmp$Assim(a, ee) / (cmp.const * a)
expect_that(sapply(ee, p$assimilation_leaf),
            equals(cmp.assimilation.leaf))

## Now, integrate to get whole plant assimilation over the light
## environment and leaf distribution.  This is probably overkill from
## a testing point of view, but this is the most complicated thing
## that happens.
cmp.assimilation.plant <- cmp$assimilation.plant(h, light.env)

f1 <- Vectorize(function(x)
                a * p$assimilation_leaf(light.env(x)) * p$q(x))
f2 <- Vectorize(function(x)
                a * p$assimilation_leaf(light.env(p$Qp(x))))
expect_that(integrate(f1, 0, h)$value,
            equals(cmp.assimilation.plant / cmp.const))
expect_that(integrate(f2, 0, 1)$value,
            equals(cmp.assimilation.plant / cmp.const,
                   tolerance=1e-7))

## Do this assimilation within the Plant class, rather than by hand.
env <- new(Spline)
env$init(hh, ee)

## Last time by hand, but using target distribution in the Plant class
f.p <- Vectorize(function(x) p$compute_assimilation_x(x, env))
expect_that(integrate(f.p, 0, h)$value * a, # leaf area
            equals(cmp.assimilation.plant / cmp.const))

## And then whole plant assimilation completely:
expect_that(p$compute_assimilation(env),
            equals(cmp.assimilation.plant / cmp.const))

## Now, look at all physiological variables:

## Compute the physiological variables and retreive them.
p$compute_vars_phys(env)
p.phys <- p$vars_phys

## 1. Assimilation:
expect_that(p.phys[["assimilation"]],
            equals(cmp.assimilation.plant / cmp.const))

## 2. Respiration:
cmp.respiration <- cmp$respiration.given.height(cmp$traits, h)
expect_that(p.phys[["respiration"]],
            equals(cmp.respiration / cmp.const))

## 3. Turnover:
cmp.turnover <- cmp$turnover.given.height(cmp$traits, h)
expect_that(p.phys[["turnover"]],
            equals(cmp.turnover))

## 4. Net production:
cmp.net.production <- cmp$net.production(cmp$traits, h, light.env)
expect_that(p.phys[["net_production"]],
            equals(cmp.net.production, tolerance=1e-7))

## 5. Reproduction fraction
cmp.reproduction.fraction <-
  cmp$ReproductiveAllocation(cmp$traits$hmat,h)
expect_that(p.phys[["reproduction_fraction"]],
            equals(cmp.reproduction.fraction))

## 6. Fecundity rate
cmp.fecundity.rate <- cmp$fecundity.rate(cmp$traits, h, light.env)
expect_that(p.phys[["fecundity_rate"]],
            equals(cmp.fecundity.rate, tolerance=1e-7))

## 7. Fraction of whole plant (mass) growth that is leaf.
cmp.leaf.fraction <- cmp$leaf.fraction(cmp$traits, h)
expect_that(p.phys[["leaf_fraction"]],
            equals(cmp.leaf.fraction))

## 8. Growth rate for leaf mass
cmp.growth.rate <- cmp$growth.rate(cmp$traits, h, light.env)
expect_that(p.phys[["growth_rate"]],
            equals(cmp.growth.rate, tolerance=1e-7))

## 9. Mortality rate
cmp.mortality.rate <- cmp$mortality.rate(cmp$traits, h, light.env)
expect_that(p.phys[["mortality_rate"]],
            equals(cmp.mortality.rate))

## Seed stuff:
seed <- new(Plant, s)

## Check that our root-finding succeeded and the leaf mass is correct:
expect_that(seed$vars_size[["mass_total"]],
            equals(pars.s$s, tolerance=1e-7))

## Check that the height at birth is correct (this happens to be how
## the R version computes size).  This is actually quite inaccurate,
## but I think that's just some instability creeping in from the
## different ways that the values involved are computed.
expect_that(seed$vars_size[["height"]],
            equals(cmp$height.at.birth(cmp$traits), tolerance=1e-4))
expect_that(seed$vars_size[["mass_leaf"]],
            equals(cmp$leaf.mass.at.birth(cmp$traits), tolerance=1e-4))

## Then, check the germination probabilities in the current light
## environment:
expect_that(seed$germination_probability(env),
            equals(cmp$germination.probability(cmp$traits, light.env),
                   tolerance=1e-5))

## Grow the plant in a constant environment
derivs <- function(t, y, pars) {
  plant <- pars$plant
  light.env <- pars$light.env
  
  plant$ode_values_set(y)
  plant$compute_vars_phys(light.env)
  plant$ode_rates
}

## Make a bigger light environment, so the plant has room to grow
env2 <- new(Spline)
hh2 <- seq(0, pars.s$hmat * 1.2, length=300)
yy2 <- light.env(hh2)
env2$init(hh2, yy2)

## Check the derivative calculations are correct
y <- c(pi, 0, 0)
pars.derivs <- list(plant=p, light.env=env2)
tmp <- derivs(t, y, pars.derivs)
p.derivs <- p.phys[c("growth_rate", "mortality_rate", "fecundity_rate")]
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
            equals(cmp.run))

## Delete the plant -- should not crash.
rm(p)
gc()

## Test copy and assignment on standalone and non-standalone plants
expect_that(tree_module$test_plant(s, FALSE, FALSE), is_true())
expect_that(tree_module$test_plant(s, TRUE,  FALSE), is_true())
expect_that(tree_module$test_plant(s, FALSE, TRUE),  is_true())
expect_that(tree_module$test_plant(s, TRUE,  TRUE),  is_true())
