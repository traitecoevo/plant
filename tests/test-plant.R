source("helper-tree.R")

context("Plant")

cmp <- make.falster()
pars.cmp <- cmp$get_parameters()

s <- new(Strategy)
pars.s <- s$get_parameters()

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

## Delete the plant -- should not crash.
rm(p)
gc()
