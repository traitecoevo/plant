## The stochastic runner for p, with each species' arrivals drawn (in C++,
## from the strategy's birth rate and ctrl$stochastic_patch_area) and `events` (an Events,
## or NULL) beside them. As in the deterministic solver, an events list naming
## introductions replaces the drawn arrivals: one arrival per time, and none
## drawn.
stochastic_runner <- function(p, env, ctrl, events = NULL) {
  ## A p saved before the area moved to Control carries patch_area. The C++
  ## converter ignores it at 1, the deterministic value, but 1 was also the old
  ## stochastic default: ignored here, such a run would draw a hundred times the
  ## plants without a word.
  if ("patch_area" %in% names(p)) {
    stop("`p$patch_area` is retired. Set a stochastic patch's area with ",
         "control(stochastic_patch_area = ) (100 m2 by default; 1 repeats ",
         "the old default), and remove patch_area", call. = FALSE)
  }
  types <- extract_RcppR6_template_types(p, "Parameters")
  obj <- do.call('StochasticPatchRunner', types)(p, env, ctrl)
  given <- !is.null(events) && length(events$time) > 0L
  if (given && any(events$type == "node_introduction")) {
    sched <- Schedule(length(p$strategies))
    sched$max_time <- p$max_patch_lifetime
    sched$add_events(events)
    obj$schedule <- sched
    return(obj)
  }
  if (given) {
    sched <- obj$schedule
    sched$add_events(events)
    obj$schedule <- sched
  }
  obj$draw_arrivals()
  obj
}

## Runs `obj` to the end, keeping the patch at its start and after every step: a
## list(time, species) where each species is a state matrix carrying an
## `is_alive` attribute. Read from obj$patch$state (the runner has no `state`
## accessor of its own; see issue #498).
collect_stochastic_steps <- function(obj) {
  res <- list(obj$patch$state)
  while (!obj$complete) {
    obj$run_next()
    res <- c(res, list(obj$patch$state))
  }
  res
}

## The collected steps as arrays: time, species and env.
tidy_stochastic_steps <- function(res) {
  time <- sapply(res, "[[", "time")
  ## `env` is what StochasticPatch reports it as, and what a deterministic
  ## run's tidy() calls it. This read used to ask for "light_env", a name
  ## nothing produced, so every element came back NULL.
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

  # offspring_production is not defined for the finite-population model (no
  # patch-density integral); omitted. patch_density likewise (see note above).
  list(time = time, species = species, env = env)
}
