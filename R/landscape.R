fitness_landscape <- function(trait_matrix, p, raw_seed_rain=FALSE) {
  n_residents <- length(p$strategies)
  if (n_residents == 0L &&
      length(p$cohort_schedule_ode_times) == 0L) {
    p$cohort_schedule_ode_times <-
      c(p$cohort_schedule_times_default, p$cohort_schedule_max_time)
  }

  p_with_mutants <- expand_parameters(trait_matrix, p)

  ebt <- run_ebt(p_with_mutants,
                 use_ode_times=length(p$cohort_schedule_ode_times) > 0)
  seed_rain <- ebt$seed_rains
  if (n_residents > 0L) {
    seed_rain <- seed_rain[-seq_len(n_residents)]
  }
  if (raw_seed_rain) {
    seed_rain
  } else {
    log(seed_rain)
  }
}

expand_parameters <- function(trait_matrix, p, mutant=TRUE) {
  if (length(mutant) != 1L) {
    stop("mutant must be scalar")
  }
  extra <- strategy_list(trait_matrix, strategy=p$strategy_default)
  n_extra <- length(extra)

  ret <- p <- validate(p) # Ensure times are set up correctly.
  ret$strategies <- c(p$strategies, extra)
  ret$is_resident <- c(p$is_resident, rep(!mutant, n_extra))
  ret$seed_rain <- c(p$seed_rain, rep(1.0, n_extra))

  ## Introduce mutants at all unique times:
  if (length(p$strategies) == 0L || !mutant) {
    times_new <- p$cohort_schedule_times_default
  } else {
    times_new <- unique(sort(unlist(p$cohort_schedule_times)))
  }
  ret$cohort_schedule_times <- c(p$cohort_schedule_times,
                                 rep(list(times_new), n_extra))

  ## Clear this if it's present:
  attr(ret, "seed_rain_out") <- NULL

  ret
}

remove_residents <- function(p) {
  if (length(p$strategies) > 0L) {
    p$strategies <- list()
    p$is_resident <- logical(0)
    p$seed_rain <- numeric(0)
    p$cohort_schedule_times <- list()
  }
  p
}
