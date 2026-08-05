## This might yet move into the compiled code if we need it there.
##
## Generate a vector of arrival times.
##
## This will be slow, but fairly easy to get right.
##' @importFrom stats rexp rpois runif splinefun
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
    ## patch_area must be passed by name: it is the *fourth* argument of
    ## stochastic_arrival_times(), and positionally it would land in `delta_t`.
    ## That silently dropped area-scaling of the arrival rate (the expected
    ## number of arrivals came out as max_time * birth_rate, whatever the area)
    ## and set the binning interval to the area, which for a large patch left
    ## only a handful of bins for a variable birth rate to be averaged over.
    times <- stochastic_arrival_times(max_time, species, patch_area = patch_area)
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
run_stochastic_collect <- function(p, env = NULL, 
                                   ctrl = Control(), 
                                   random_schedule=TRUE) {
  
  types <- extract_RcppR6_template_types(p, "Parameters")

  if (is.null(env)) {
    env <- Environment(types[[1]])
  }
  
  ## Per-step snapshot of the patch: a list(time, species) where each species is
  ## a state matrix carrying an `is_alive` attribute. Read from obj$patch$state
  ## (the runner has no `state` accessor of its own; see issue #498).
  collect <- function(obj) {
    obj$patch$state
  }
  types <- extract_RcppR6_template_types(p, "Parameters")
  obj <- do.call('StochasticPatchRunner', types)(p, env, ctrl)
  if (random_schedule) {
    obj$node_schedule <- stochastic_schedule(p)
  }

  res <- list(collect(obj))

  while (!obj$complete) {
    obj$run_next()
    res <- c(res, list(collect(obj)))
  }

  time <- sapply(res, "[[", "time")
  ## `env` is what StochasticPatch reports it as, and what run_scm's collected
  ## output calls it. This read used to ask for "light_env", a name nothing
  ## produced, so every element came back NULL.
  env <- lapply(res, "[[", "env")
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
              env=env,
              # offspring_production is not defined for the finite-population
              # model (no patch-density integral); omitted. patch_density
              # likewise (see note above).
              p=p)

  ret
}
