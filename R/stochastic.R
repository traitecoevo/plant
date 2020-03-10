## This might yet move into the compiled code if we need it there.
##
## Generate a vector of arrival times.
##
## This will be slow, but fairly easy to get right.
##' @importFrom stats rexp
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
  for (i in seq_along(seed_rain)) {
    sched$set_times(stochastic_arrival_times(max_time, seed_rain[[i]]), i)
  }
  sched
}

##' Run a stochastic simulation of a patch, given a Parameters
##'
##' This one might need to be made differently so that different
##' schedules can be added easily.
##' Not sure if this is how we will generally want to do this.
##' Consider this function liable to change.
##'
##' @title Run a stochastic patch, Collecting Output
##' @param p A \code{\link{FF16_Parameters}} object
##' @param random_schedule setting to TRUE causes algorithm to generate
##' a random schedule based on seed rain and area.
##' @author Rich FitzJohn
##' @export
run_stochastic_collect <- function(p, random_schedule=TRUE) {
  collect <- function(obj) {
    obj$state
  }
  types <- extract_RcppR6_template_types(p, "Parameters")
  obj <- do.call('StochasticPatchRunner', types)(p)
  if (random_schedule) {
    obj$schedule <- stochastic_schedule(p)
  }

  res <- list(collect(obj))

  while (!obj$complete) {
    obj$run_next()
    res <- c(res, list(collect(obj)))
  }

  time <- sapply(res, "[[", "time")
  light_env <- lapply(res, "[[", "light_env")
  species <- lapply(res, "[[", "species")

  ## The aperm() here means that dimensions are
  ## [variable,time,plant], so that taking species[[1]]["height",,]
  ## gives a matrix that has time down rows and plants across columns
  ## (so is therefore plottable with matplot)

  n_spp <- length(species[[1]])

  species_is_alive <- lapply(seq_len(n_spp), function(i)
    t(pad_list_to_array(lapply(species, function(x) attr(x[[i]], "is_alive")))))
  species <- lapply(seq_len(n_spp), function(i)
    aperm(pad_list_to_array(lapply(species, "[[", i)), c(1, 3, 2)))
  attr(species, "is_alive") <- species_is_alive

  patch_density <- obj$patch$environment$disturbance_regime$density(time)

  ret <- list(time=time,
              species=species,
              light_env=light_env,
              seed_rain=obj$seed_rains,
              patch_density=patch_density,
              p=p)

  ret
}
