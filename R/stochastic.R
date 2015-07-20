## This might yet move into the compiled code if we need it there.
##
## Generate a vector of arrival times.
##
## This will be slow, but fairly easy to get right.
stochastic_arrival_times <- function(max_time, seed_rain_total, n=NULL) {
  if (is.null(n)) {
    n <- max_time * seed_rain_total
  }
  ret <- numeric(0)
  t0 <- 0.0
  while (t0 < max_time) {
    t <- t0 + cumsum(rexp(n, seed_rain_total))
    t0 <- t[[n]]
    if (t0 > max_time) {
      t <- t[t <= max_time]
    }
    ret <- c(ret, t)
  }
  ret
}

stochastic_schedule <- function(p) {
  max_time  <- p$cohort_schedule_max_time
  seed_rain <- p$seed_rain * p$patch_area
  n_species <- length(seed_rain)
  sched <- CohortSchedule(n_species)
  sched$max_time <- max_time
  for (i in seq_len(seed_rain)) {
    sched$set_times(stochastic_arrival_times(max_time, seed_rain[[i]]), i)
  }
  sched
}
