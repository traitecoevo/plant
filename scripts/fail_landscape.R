
# Documents case when landscape function does not give correct fitness

library(tree)

time.disturbance <- 5
strategy <- new(Strategy, list(c_r1=0.5,c_r2=0))

parameters <- new(Parameters)
parameters$set_parameters(list(patch_area=1.0))
parameters$set_control_parameters(fast.control())
parameters$disturbance <- new(Disturbance, time.disturbance)

max.t <- parameters$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(max.t)

sys <- list(traits=cbind(lma=0.1978791),
            seed_rain=505.55407,
            times=list(times0))

build.args <- list(nsteps=10, eps=1e-3, verbose=TRUE)


## update population dynamics
traits <- sys[["traits"]]
p <- parameters$copy()
for (i in seq_len(nrow(traits))){
  new.strategy <- strategy$copy()
  new.strategy$set_parameters(as.list(traits[i,]))
  p$add_strategy(new.strategy)
}
p$seed_rain <- sys[["seed_rain"]]

  schedule <- new(CohortSchedule, length(sys[["times"]]))
  schedule$max_time <- max.t
  schedule$all_times <- sys[["times"]]

res <- build.schedule(p, schedule, build.args$nsteps, build.args$eps,
                      progress=FALSE, verbose=build.args$verbose)

rain.out <- unname(attr(res, "seed_rain", exact=TRUE)[,"out"])

# Here's the ratio of seed rain in to out
rain.out/ sys[["seed_rain"]]

# Should be similar to fitness estimated for resident but isn't

schedule$ode_times <- attr(res, "ebt")$ode_times
landscape("lma",  sys[["traits"]][,"lma"], p, schedule)
