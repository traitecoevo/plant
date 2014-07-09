# Documents case when landscape function does not give correct fitness

library(tree)

time.disturbance <- 5
strategy <- new(Strategy, list(c_r1=0.5, c_r2=0))

parameters <- new(Parameters)
parameters$set_parameters(list(patch_area=1.0))
parameters$set_control_parameters(fast.control())
parameters$disturbance <- new(Disturbance, time.disturbance)
parameters$set_control_parameters(list(schedule_verbose=TRUE))

sys <- list(traits=cbind(lma=0.1978791),
            seed_rain=505.55407)

## update population dynamics
traits <- sys[["traits"]]
p <- parameters$copy()
for (i in seq_len(nrow(traits))){
  new.strategy <- strategy$copy()
  new.strategy$set_parameters(as.list(traits[i,]))
  p$add_strategy(new.strategy)
}
p$seed_rain <- sys[["seed_rain"]]

res <- build_schedule(p)

rain_out <- unname(attr(res, "seed_rain", exact=TRUE)[,"out"])

# Here's the ratio of seed rain in to out
rain_out / sys[["seed_rain"]]

# Should be similar to fitness estimated for resident but isn't
(w1 <- landscape("lma",  sys[["traits"]][,"lma"], p, res))

# works when we pass in strategy
p$strategy_default <- strategy
(w2 <- landscape("lma",  sys[["traits"]][,"lma"], p, res))
