library(tree)

## Copied from direct-tree-fun.R
add_strategy_trait <- function(trait, value, p, s, seed_rain=NULL) {
  if (p$size > 0) {
    stop("Must have empty environment")
  }
  p <- p$copy()
  s <- s$copy()
  s$set_parameters(structure(list(value), names=trait))
  p$add_strategy(s)
  if (!is.null(seed_rain)) {
    p$seed_rain <- seed_rain
  }
  p
}

add_ode_times <- function(trait, value, p, schedule, s, seed_rain=1) {
  p <- add_strategy_trait(trait, value, p, s, seed_rain)
  ebt <- new(EBT, p)
  ebt$cohort_schedule <- schedule$copy()
  ebt$run()
  ebt$ode_times
}

p <- new(Parameters)
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)
schedule <- schedule.from.times(times0)

hmat.r <- c(1, 50) # only fails when lower bound is 1, not 2.
trait <- "hmat"
values <- seq_log(hmat.r[[1]], hmat.r[[2]], 20)

s <- new(Strategy)
if (length(schedule$ode_times) == 0) {
  values.mid <- exp(mean(log(range(values))))
  schedule$ode_times <- add_ode_times(trait, values.mid, p, schedule, s)
}

values <- 1
p.empty <- p$copy()
p.empty$clear()
for (i in values) {
  l <- structure(list(i), names=trait)
  p.empty$add_strategy_mutant(new(Strategy, l))
}

schedule.empty <- new(CohortSchedule, length(values))
schedule.empty$max_time  <- schedule$max_time
schedule.empty$ode_times <- schedule$ode_times
schedule.empty$all_times <-
  rep(list(unique(sort(unlist(schedule$all_times)))), length(values))
## ebt.empty <- run.ebt(p.empty, schedule.empty)

matrix(ebt$ode_values, 4)[1,]

ebt <- new(EBT, p.empty)
ebt$cohort_schedule <- schedule.empty
ebt$run()

patch <- ebt$patch
matrix(patch[[5]]$ode_values, 4)[1,]
matrix(patch[[4]]$ode_values, 4)[1,]

m <- matrix(ebt$ode_values, 4)
rownames(m) <- c("height", "mortality", "fecundity", "log_density")
r <- matrix(ebt$ode_rates, 4)
rownames(r) <- c("height", "mortality", "fecundity", "log_density")

ebt <- new(EBT, p.empty)
ebt$cohort_schedule <- schedule.empty
for (


while (!ebt$complete) {
  ebt$run_next()
}


sched <- ebt$cohort_schedule
sched$remaining
