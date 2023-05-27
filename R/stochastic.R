## This might yet move into the compiled code if we need it there.
##
## Generate a vector of arrival times.
##
## This will be slow, but fairly easy to get right.
##' @importFrom stats rexp rpois runif
stochastic_arrival_times <- function(max_time, species, delta_t = 0.1, patch_area = 1) {
  ret <- numeric(0)
  t0 <- 0.0
  t1 <- t0 + delta_t

  # first calculate average arrival rate in this interval
  if (species$is_variable_birth_rate) {
    interpolated <- splinefun(species$birth_rate_x, species$birth_rate_y)
  } else {
    rate <- species$birth_rate_y
  }

  while (t0 < max_time) {
    # first calculate average arrival rate in this interval
    if (species$is_variable_birth_rate) {
      x = seq(t0, t1, len = 10)
      rate = mean(interpolated(x))
    }

    # now calculate actual number arriving, given the rate
    n <- rpois(1,  delta_t * rate * patch_area)

    # now generate arrival times in this interval, from a uniform distribution
    t <- sort(runif(n, t0, t1))
    ret <- c(ret, t)

    # update for next iteration
    t0 <- t1
    t1 <- t1 + delta_t
  }

  ret[ret < max_time]
}


stochastic_schedule <- function(p) {
  patch_area <- p$patch_area
  max_time  <- p$max_patch_lifetime
  n_species <- length(p$strategies)

  sched <- NodeSchedule(n_species)
  sched$max_time <- max_time

  for (i in 1:n_species) {
    species <- p$strategies[[i]]
    times <- stochastic_arrival_times(max_time, species, patch_area)
    sched$set_times(times, i)
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
##' @param env Environment object
##' @param ctrl Control object
##' @param random_schedule setting to TRUE causes algorithm to generate
##' a random schedule based on offspring arrival and area.
##' @author Rich FitzJohn
##' @export
run_stochastic_collect <- function(p, env = make_environment(parameters = p), 
                                   ctrl = scm_base_control(), 
                                   random_schedule=TRUE) {
  collect <- function(obj) {
    obj$state
  }
  types <- extract_RcppR6_template_types(p, "Parameters")
  obj <- do.call('StochasticPatchRunner', types)(p, env, ctrl)
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

  # Not sure how this works in a stochastic patch, I wouldn't have thought
  # meta-populations worked without a continuous gradient to integrate over.
  # patch_density <- obj$patch$density(time)

  ret <- list(time=time,
              species=species,
              light_env=light_env,
              offspring_production=obj$offspring_production,
              # patch_density=patch_density,
              p=p)

  ret
}
