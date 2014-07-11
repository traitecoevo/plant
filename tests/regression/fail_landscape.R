## Documents case when landscape function did not give correct
## fitness, because the strategy that landscape built off was not the
## same as the resident.
library(tree)

time_disturbance <- 5
strategy <- new(Strategy, list(c_r1=0.5, c_r2=0))

parameters <- new(Parameters)
parameters$set_parameters(list(patch_area=1.0))
parameters$set_control_parameters(fast_control())
parameters$disturbance <- new(Disturbance, time_disturbance)
parameters$set_control_parameters(list(schedule_verbose=TRUE))

sys <- list(traits=cbind(lma=0.1978791),
            seed_rain=505.55407)

## Set up resident population.  The problem here is that this is based
## on `strategy`, not on the strategy in
## `parameters$strategy_default`.  If we'd used
## `expand_parameters(..., mutant=FALSE)`, the base strategy would be
## `$strategy_default`, and we'd not have this issue.
traits <- sys[["traits"]]
resident <- strategy$copy()
resident$set_parameters(as.list(traits[,1]))
p <- parameters$copy()
p$add_strategy(resident)
p$seed_rain <- sys[["seed_rain"]]

res <- build_schedule(p)

rain_out <- unname(attr(res, "seed_rain", exact=TRUE)[,"out"])

## Here's the ratio of seed rain in to out: should be fairly close to
## 1:
rain_out / sys[["seed_rain"]]

## Should be similar to fitness estimated for resident but doesn't
## appear to be, because the mutant is constructed with a different
## base strategy.
(w1 <- landscape("lma",  sys[["traits"]][,"lma"], p, res))

## Things work as exected when we set the default strategy though.
p$strategy_default <- strategy
(w2 <- landscape("lma",  sys[["traits"]][,"lma"], p, res))

library(testthat)
expect_that(w2, equals(rain_out / sys[["seed_rain"]],
                       tolerance=1e-5))
