source("helper-tree.R")

## TODO: Write a wrapper function within the Falster R model that
## collects all the physiological bits together.

context("Plant")

cmp <- make.falster()
pars.cmp <- cmp$get_parameters()

s <- new(Strategy)
core <- list(lma=0.1978791, hmat=16.5958691, rho=608, s=3.8e-5)
s$set_parameters(core)
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

yy.C <- sapply(ee, p$assimilation_leaf)
yy.R <- cmp$Assim(a, ee) / (cmp.const * a)

expect_that(yy.C, equals(yy.R))

## Now, integrate, using both the C++ and R versions:
f.C <- Vectorize(function(x)
                 p$assimilation_leaf(light.env(x)) * p$q(x))
f.R <- function(x) cmp$leaf.pdf(x, h) * cmp$Assim(a, light.env(x)) /
  (cmp.const* a)

Y.C <- integrate(f.C, 0, h)
Y.R <- integrate(f.R, 0, h)
expect_that(Y.C, equals(Y.C))

## Also can do this with the alternative approach:
g.C <- Vectorize(function(x) p$assimilation_leaf(light.env(p$Qp(x))))
g.R <- function(x) cmp$Assim(a, light.env(cmp$leaf.icdf(x, h))) /
  (cmp.const * a)

expect_that(integrate(g.C, 0, 1)$value,
            equals(Y.C$value, tolerance=1e-7))
expect_that(integrate(g.R, 0, 1)$value,
            equals(Y.R$value, tolerance=1e-7))

#### Messy below here.

env <- new(Spline)
env$init(hh, ee)
tmp <-
  integrate(Vectorize(function(x) p$compute_assimilation_x(x, env)),
            0, h)$value
expect_that(tmp, equals(Y.C$value))

## TODO: Need to actually work out how to compute light averaged
## photosynthesis in Daniel's R version.  To be consistent it has a
## different constant...
cmp.assimilation <- Y.C$value * a * cmp.const

expect_that(p$compute_assimilation(env),
            equals(cmp.assimilation / cmp.const))

p$compute_vars_phys(env)
## This shows obvious issues; respiration far too high (seen that
## before).  turover and net_production failed to calculate.  Others
## are understandible.
ans <- p$vars_phys

expect_that(ans[["assimilation"]],
            equals(cmp.assimilation / cmp.const))

cmp.respiration <- 
  cmp$Respiration(size.p[["leaf_area"]],
                  size.p[["mass_sapwood"]]/cmp$traits$rho,
                  size.p[["mass_bark"]]/cmp$traits$rho,
                  size.p[["mass_root"]])
expect_that(ans[["respiration"]],
            equals(cmp.respiration / cmp.const))

cmp.turnover <- 
  cmp$Turnover(cmp$traits,
               size.p[["mass_leaf"]],
               size.p[["mass_sapwood"]],
               size.p[["mass_bark"]],
               size.p[["mass_root"]])
expect_that(ans[["turnover"]],
            equals(cmp.turnover))

## Like the Assimilation one, we actually can't do this directly as
## the light environment does not integrate?
## TODO: where does the difference here come from?  Higher order than
## I would have thought.  Looks like it's from the integrator.  Check
## the error bounds more carefully.
cmp.net.production <-
  cmp.assimilation - cmp.respiration - cmp.turnover
expect_that(ans[["net_production"]],
            equals(cmp.net.production, tolerance=1e-7))

## Actually identical.
cmp.reproduction.fraction <- cmp$ReproductiveAllocation(cmp$traits$hmat,h)
expect_that(ans[["reproduction_fraction"]],
            equals(cmp.reproduction.fraction))

## TODO: Reproduction not computed in the R version, apparently, just
## the fraction.  For now here is a constant.
expect_that(ans[["fecundity_rate"]],
            equals(0.530234597058596))

## TODO: Not actually computed in R version.  Can get dMtdt and things
## like that but not exactly.
## 
## TODO: fraction is really bad as it should be 0-1 if it really is
## one.
##
## TODO: given that cmp.net.production * (1-cmp.reproduction.fraction)
## is 10, I don't see how this is 15.
## expect_that(ans[["leaf_fraction"]],
##             equals())

## TODO: Can't actually do this directly because Production does not
## use environment.  Also, the dmldmt is not computed.
## cmp.growth.rate <-
##   cmp.net.production * (1-cmp.reproduction.fraction)

 
## expect_that(ans[["growth_rate"]],
##             equals(cmp.growth.rate))



## Delete the plant -- should not crash.
## rm(p)
## gc()

