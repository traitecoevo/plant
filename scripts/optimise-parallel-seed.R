devtools::load_all()

library(doRedis)
setPackages(c("tidyverse", "laGP", "lhs"))

source("scripts/optimise-bayesopt.R")



# calibrate against tree yield formula for a given max. bio
tyf <- function(t, m = 50, g=12.25, r=1){
  k = 2 * g - 1.25
  agb = r * m * exp(-k / t)
  # bgb = agb * .67
  
  return(agb / 10) # t.ha-1 -> kg.m-2
}

setExport(ls())


# parameter bounds
bounds = list(B_lf1 = list(min = 0.5, max = 1),
              hmat = list(min = 4, max = 8),
              eta = list(min = 5, max = 12),
              k_2 = list(min = 0, max = 0.5),
              birth_rate = list(min = 0.01, max = 2))



# generate seed
p <- create_schedule(max_patch_lifetime = 50, optimise_schedule = T)

# register queue
registerDoRedis('jobs')

# add workers
n_threads = 16
startLocalWorkers(n=n_threads, queue="jobs", linger=1)

# load plant
foreach(i = 1:n_threads) %dopar%
    devtools::load_all()

    
# dispatch work - will see how 'job::job' goes for detatched runs
# also keen to use external cache for results
setProgress(TRUE)
seed <- generate_seeds_parallel(calibrate, bounds, n = 100, target = tyf, parameters = p)

#
removeQueue('jobs')

