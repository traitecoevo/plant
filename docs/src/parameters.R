## ---
## title: "plant: A package for modelling forest trait ecology & evolution: _Modifying parameters_"
## ---

library(plant)

## There are a large number of parameters to the physiological model,
## but all are changeable.  The default strategy is detailed in the
## "physiology" vignette:
s <- FF16_Strategy()
names(s)

## The strategy object here is a special object of class `r class(s)`.
## In contrast with most of the "reference" objects used by `plant`,
## this is simply a list with a class attribute.  Some validation will
## be done on the parameters every time it is passed into the C++ bits
## of the model.

## All parameters except for `control` are floating point (decimal)
## numbers:
unlist(s[names(s) != "control"])

## All of these parameters can be directly changed.  For example, we
## can double the leaf mass per unit leaf area (lma) value:
s$lma <- s$lma * 2

## and then from this construct a plant:
pl <- FF16_Plant(s)
pl$strategy$lma

## `lma` affects a few places in the model; see [the source
## code](https://github.com/traitecoevo/plant/blob/master/src/ff16_strategy.cpp),
## but only really for converting from leaf area to leaf mass (leaf
## mass being leaf area multiplied by leaf mass per leaf area).
## However, as a component of the leaf economic spectrum we imagine
## LMA as affecting a number of other components of the model.

## We capture this through what we call a "hyper-parameterisation";
## additional parameters and functions that mean that changing one
## parameter might affect a number of other lower-level parameters.
## Our default hyper-parameterisation is via the
## `FF16_hyperpar` function:
FF16_hyperpar

## You will rarely need to call this function directly (see below) but
## note how setting of lma affects parameters `k_l`, `a_p1` and `r_l`
FF16_hyperpar(trait_matrix(0.1, "lma"), s)

## These are:
##   * `k_l`: Turnover rate for leaves
##   * `a_p1`: Leaf photosynthesis per area
##   * `r_l`: Leaf respiration per mass

## This means that it is possible to implement different trade-offs
## between parameters relatively easily by modifying the
## hyper-parameterisation function in R, rather than having to modify
## the underlying physiological model in C++.

## The `trait_matrix` function is a simple wrapper that just makes
## sure the trait matrix has the right format:
trait_matrix(0.1, "lma")
trait_matrix(c(0.1, 500), "rho")

## To make use of the hyper-parameterisation, the preferred way of
## setting parameters is through the utility functions `strategy` and
## `strategy_list`.  These take a `Parameters` object:
p <- FF16_Parameters()

## The Parameters object mostly contains information about the patch:
names(p)

## and it comes pre-set with the hyper-parameterisation function:
identical(p$hyperpar, FF16_hyperpar)

## `k_I` is the light extinction coefficient,
## `disturbance_mean_interval` is the mean disturbance interval.  The
## `strategy_default` member is a Strategy:
class(p$strategy_default)

## This is the Strategy object that all others will be built from by
## difference.  Running
s <- strategy(trait_matrix(0.1, "lma"), p)

## will create a strategy `s` where lma is set but also all the
## parameters that *depend* on lma.

## The function `strategy_list` can be used to create a list of
## strategies:
lma <- trait_matrix(seq(0.1, 0.5, length.out=5), "lma")
FF16_hyperpar(trait_matrix(lma, "lma"), s)

ss <- strategy_list(lma, p)
length(ss)

## We can then use standard R commands to extract variable from this list
sapply(ss, function(x) x$lma)
sapply(ss, function(x) x$k_l)
sapply(ss, function(x) x$a_p1)
sapply(ss, function(x) x$r_l)

## There's a convenience function `plant_list` that returns a set of
## *plants* based on a vector of traits:
pp <- plant_list(lma, p)

sapply(pp, function(p) p$area_leaf_above(0))

## In addition to the physiological parameters there are large number
## of "control" parameters that affect the behaviour of the various
## numerical algorithms used (note that in contrast to the
## physiological parameters these have a variety of types)
p$control

## The defaults are rather too slow for many uses, so
## `scm_base_parameters` provides a faster set by using `fast_control`
## to set many of these to less accurate values.

p2 <- scm_base_parameters()
p2$control[unlist(p$control) != unlist(p2$control)]
