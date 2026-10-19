##' Construct a \code{Control} object. \code{control()} is a lowercase alias for
##' the \code{Control()} constructor, whose defaults are the pragmatic,
##' fast-ish settings used for essentially all of plant's runs (see
##' \code{control.cpp}). \code{control_accurate()} tightens the ODE and schedule
##' tolerances for high-accuracy runs at the cost of speed.
##'
##' @title Control presets
##' @param ... Named control fields, passed to \code{Control()}.
##' @param base An optional \code{Control} object to tighten; defaults are used
##'   if omitted.
##' @return A \code{Control} object.
##' @rdname control_presets
##' @export
control <- function(...) Control(...)

##' @rdname control_presets
##' @export
control_accurate <- function(base = Control()) {
  base$ode_tol_rel       <- 1e-6
  base$ode_tol_abs       <- 1e-6
  base$ode_step_size_max <- 1e-1
  base$schedule_eps      <- 1e-3
  base
}


##' Basic default settings for a given strategy, environment only really
##' used for templating initially and will be overloaded later by passing
##' an environment to the SCM API (suggesting perhaps the template could be
##' removed).
##' @title Basic default parameters for a given strategy
##' @author Rich FitzJohn
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @param env And environment object
##' @export
scm_base_parameters <- function(type = NA, env = environment_type(type)) {
  Parameters(type, env)()
}


##' Run the SCM.
##'
##' The node-introduction schedule can be adaptively refined in C++ by setting
##' \code{refine_schedule = TRUE} (this replaces the former \code{build_schedule}
##' function). Setting \code{collect = TRUE} returns tidied output collected at
##' every ODE step (replacing the former \code{run_scm_collect}); otherwise the
##' \code{SCM} object itself is returned for interrogation.
##'
##' @title Run SCM
##' @param p Parameters object
##' @param env Environment object (defaults to the strategy's environment)
##' @param ctrl Control object
##' @param events An \code{\link{events}} object giving the discrete events to
##'   apply during the run — rainfall pulses, harvest, and node introductions
##'   themselves. When \code{NULL} (the default) the schedule is taken from
##'   \code{p$node_schedule_times}, i.e. introductions only.
##' @param refine_schedule Should the node-introduction schedule be adaptively
##'   refined before/while running (using \code{schedule_eps} and
##'   \code{schedule_nsteps} from \code{ctrl})? Refinement records the ODE
##'   schedule its final run took into \code{p$ode_times},
##'   \code{p$ode_step_sizes}, \code{p$ode_soil_input_slopes} and
##'   \code{p$ode_soil_substep_ends}, so a later run of those parameters replays
##'   it exactly rather than choosing its own steps again.
##' @param record_trajectory Should the run keep the state at every accepted
##'   step? A gradient sweeps those states and cannot recover them from a
##'   finished run, so a run that did not keep them is repeated -- one whole
##'   forward integration. Asked for here because the flag has to be set before
##'   the run that fills it, and this function is where that run happens. With
##'   \code{refine_schedule = TRUE} the refinement's own runs do not keep them;
##'   only the final run does.
##' @param collect Should tidied results be collected at every step and
##'   returned (instead of the \code{SCM} object)?
##' @param run Run the SCM? \code{FALSE} returns it built and unrun, for a
##'   caller that drives it leg by leg: \code{\link{run_until_converged}} and
##'   \code{SCM$run_to()}. It cannot be combined with \code{collect} or
##'   \code{refine_schedule}, which are about the run.
##' @return When \code{collect = FALSE}, an \code{SCM} object. When
##'   \code{collect = TRUE}, a list of tidied patch output with
##'   \code{offspring_production}, \code{net_reproduction_ratios} and the
##'   (possibly refined) parameters \code{p}.
##' @author Rich FitzJohn
##' @rdname run_scm
##' @export
run_scm <- function(p, env = NULL,
                    ctrl = control(),
                    refine_schedule = FALSE, collect = FALSE,
                    record_trajectory = FALSE, events = NULL, run = TRUE) {
  if (!run && (collect || refine_schedule)) {
    stop("`run = FALSE` returns an unrun SCM, so `collect` and ",
         "`refine_schedule`, which are about the run, cannot be asked for")
  }

  types <- extract_RcppR6_template_types(p, "Parameters")

  if (is.null(env))
    env <- Environment(types[[1]])

  # An ODE schedule carried by the parameters is taken, because a schedule is
  # there to be used: p$ode_times with p$ode_step_sizes and the soil's records
  # replays a recorded run exactly, and p$ode_times alone stops at a grid the
  # caller chose. To integrate freely, carry neither.

  ## No events supplied: the schedule comes from p$node_schedule_times, as it
  ## did before events existed. An empty Events object is how that is signalled
  ## across the boundary.
  if (is.null(events))
    events <- empty_events()

  scm <- do.call('SCM', types)(p, env, events, ctrl)
  if (!run) {
    scm$record_trajectory <- record_trajectory
    return(scm)
  }
  if (collect) {
    scm$collect <- TRUE
  }

  if (refine_schedule) {
    # Refinement's runs are bisected against and discarded, so they keep no
    # states; the run a sweep walks is one more, after the schedule settles.
    scm$refine_schedule()
    if (record_trajectory) {
      scm$record_trajectory <- TRUE
      scm$run()
    }
  } else {
    scm$record_trajectory <- record_trajectory
    scm$run()
  }

  if (!collect) {
    return(scm)
  }

  results <- lapply(scm$history, "[[", "state") |> tidy_patch()
  results[["offspring_production"]] <- scm$offspring_production
  results[["net_reproduction_ratios"]] <- scm$net_reproduction_ratios
  results[["p"]] <- scm$parameters
  ## Events are supplied separately from `p`, so `p` alone does not describe the
  ## run. Carry both the requested schedule and what was actually applied, or a
  ## collected result silently loses the whole event record -- including how much
  ## of each pulse the soil accepted and how much it shed.
  results[["events"]] <- scm$events
  results[["event_log"]] <- scm$event_log

  results
}

##' Export the full state of a patch from a (run) \code{SCM} so it can be
##' re-imported to seed a new run (see \code{\link{set_initial_state}}). The
##' exported object captures everything needed to reproduce the patch's forward
##' trajectory: every node's ODE state plus the per-node birth bookkeeping
##' (introduction time, patch-age density and survival probability at birth)
##' that is not part of the ODE state but feeds the rates and lifetime-fitness
##' integrals, the patch age, and the not-yet-introduced ("residual") portion of
##' the node-introduction schedule. A resume re-derives the two birth weights
##' from its own regime (see \code{\link{set_initial_state}}); they are exported
##' as a record of the run that made them.
##'
##' @title Export patch state from an SCM
##' @param scm An \code{SCM} object that has been run.
##' @param step Optional 1-based index into \code{scm$history} (requires the run
##'   to have been performed with \code{collect = TRUE}); when \code{NULL}
##'   (default) the SCM's current/final patch is used.
##' @return A list describing the patch state, suitable for
##'   \code{\link{set_initial_state}}: \code{time}, \code{n} (nodes per species),
##'   \code{ode_state} (flat), the per-species lists \code{node_times},
##'   \code{patch_density} and \code{pr_patch_survival}, the residual
##'   \code{node_schedule_times}, and \code{events}, the run's events still to
##'   come -- its residual introductions and every action the run had not yet
##'   applied when the state was taken (an action at the export time itself is
##'   pending in a collected snapshot and in a run that \code{SCM$run_to()}
##'   stopped there, and applied in a finished run's final patch). Events are a \code{\link{run_scm}} argument rather than a field
##'   of \code{Parameters}, so pass \code{state$events} to the resumed run's
##'   \code{events} argument; a resumed run given none re-runs the residual
##'   introductions alone and silently loses any pulse or harvest still ahead.
##'
##' With the absorbing top layer on (\code{Control$absorb_epsilon}), each
##' species' block of \code{ode_state} ends with the layer: its pinned
##' individual's states, then its density, its two fitness weights and its two
##' seed totals. A state exported without the layer has no such block, so
##' resuming it with the layer on needs that many zeros inserted after each
##' species' nodes.
##' @seealso \code{\link{set_initial_state}}, \code{\link{run_scm}}
##' @export
export_patch_state <- function(scm, step = NULL) {
  patch <- if (is.null(step)) scm$patch else scm$history[[step]]
  time <- patch$time
  species <- patch$species
  tol <- 1e-8

  ## Residual schedule: the introductions not yet made. A seeded node carries
  ## its introduction time, so the residual is every scheduled time at or after
  ## `time` that no seeded node of that species carries. A patch snapshot sits
  ## integrated *to* its next introduction time without having introduced that
  ## node yet, so the node due at `time` itself stays in the residual; a run
  ## that ended on an introduction (one scheduled at max_patch_lifetime) has
  ## made it, and a time a seeded node carries is never residual.
  node_times <- lapply(species, function(s) s$node_times)
  residual <- Map(function(tt, made) tt[tt >= time - tol & !(tt %in% made)],
                  scm$parameters$node_schedule_times, node_times)

  list(
    time = time,
    n = as.integer(vapply(species, function(s) s$size, numeric(1))),
    ode_state = patch$ode_state,
    node_times = node_times,
    patch_density = lapply(species, function(s) s$patch_densities),
    pr_patch_survival = lapply(species, function(s) s$pr_patch_survival_at_birth),
    node_schedule_times = residual,
    events = residual_events(scm, time, residual, tol,
                             final = isTRUE(scm$complete) &&
                               (is.null(step) || step == length(scm$history)))
  )
}

## The events a resumed run still has to apply: the residual introductions, and
## every action the run had not applied when this state was taken. An action
## strictly after `time` is pending. An action at exactly `time` is pending in a
## collected snapshot, which is taken on arrival at an instant before the
## instant's actions run, and in a run that run_to() stopped there, which stops
## on arrival too; it is applied in a finished run's final patch, which the run
## leaves only after applying everything scheduled at its end. The event log
## cannot decide this, because it covers the whole run whichever snapshot is
## exported.
residual_events <- function(scm, time, residual, tol, final) {
  ev <- scm$events
  n <- length(ev$time)
  if (n == 0L) {
    return(ev)
  }
  keep <- logical(n)
  for (i in seq_len(n)) {
    if (ev$type[[i]] == "node_introduction") {
      tt <- residual[[ev$target_index[[i]]]]
      keep[i] <- any(abs(tt - ev$time[[i]]) <= tol)
    } else if (ev$time[[i]] > time + tol) {
      keep[i] <- TRUE
    } else if (abs(ev$time[[i]] - time) <= tol) {
      keep[i] <- !final
    }
  }
  Events(time = ev$time[keep], type = ev$type[keep], target = ev$target[keep],
         target_index = ev$target_index[keep], params = ev$params[keep])
}

##' Build an initial patch state from a specified size distribution, for seeding
##' a patch at age 0 with pre-existing plants instead of growing it from empty
##' (the ecological motivation of \code{plant} issue #304). The returned object
##' is consumed by \code{\link{set_initial_state}}.
##'
##' Each species' initial nodes are described by their heights and densities;
##' the remaining ODE state (mortality, fecundity, accumulated reproduction,
##' heartwood, ...) starts at zero, all nodes are introduced at patch age 0, and
##' the recruitment schedule continues for \code{t > 0} (the seeded distribution
##' replaces the \code{t = 0} recruit). Pathologically large/dense initial
##' conditions can produce non-finite densities; \code{run_scm} guards against
##' this and errors with a suggestion to use more plausible inputs.
##'
##' @title Build an initial size distribution
##' @param p A \code{Parameters} object.
##' @param heights Per-species node heights: a numeric vector (single species) or
##'   a list of numeric vectors (one per species).
##' @param densities Per-species node densities (same shape as \code{heights}).
##'   Supply this or \code{log_densities}.
##' @param log_densities Per-species node log-densities (alternative to
##'   \code{densities}).
##' @param env Environment object (defaults to the strategy's environment).
##' @param ctrl Control object.
##' @param birth_dates Per-species node birth dates, before 0 and ascending as
##'   the heights descend (same shape as \code{heights}). Required on the
##'   birth-date coordinate, where each node's density is its birth rate times
##'   its survival and \code{densities} are densities in birth date: the
##'   survival is set to give that density, and each interval's establishment to
##'   its width, as if every seed established. Elsewhere every node is born at 0.
##' @return A state list suitable for \code{\link{set_initial_state}}.
##' @seealso \code{\link{set_initial_state}}, \code{\link{export_patch_state}}
##' @export
make_initial_state <- function(p, heights, densities = NULL,
                               log_densities = NULL, env = NULL,
                               ctrl = control(), birth_dates = NULL) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  if (is.null(env)) {
    env <- Environment(types[[1]])
  }
  if (ctrl$absorb_epsilon > 0) {
    stop("make_initial_state() lays out the nodes' state by hand and has no ",
         "absorbing top layer; build the state with control$absorb_epsilon = 0")
  }
  n_spp <- length(p$strategies)

  as_list <- function(z) if (is.list(z)) z else list(z)
  heights <- as_list(heights)
  if (length(heights) != n_spp) {
    stop("`heights` must have one entry per species (", n_spp, ")")
  }
  if (is.null(log_densities)) {
    if (is.null(densities)) {
      stop("supply either `densities` or `log_densities`")
    }
    log_densities <- lapply(as_list(densities), log)
  } else {
    log_densities <- as_list(log_densities)
  }

  ## Learn the node ODE layout and the (node-free) environment ODE tail from a
  ## fresh patch, plus the patch-age-0 disturbance weights for birth bookkeeping.
  patch <- do.call("Patch", types)(p, env, ctrl)
  ode_names <- patch$species[[1]]$new_node$ode_names
  hi <- match("height", ode_names)
  birth_date <- ctrl$node_density_in_birth_date
  ldi <- match(if (birth_date) "mortality" else "log_density", ode_names)
  ii <- match(c("interval_establishment", "interval_establishment_moment"),
              ode_names)
  if (is.na(hi) || is.na(ldi) || (birth_date && anyNA(ii))) {
    stop("could not locate the node's height and density in its ODE names")
  }
  if (birth_date) {
    if (is.null(birth_dates)) {
      stop("the birth-date coordinate needs each node's `birth_dates`")
    }
    birth_dates <- as_list(birth_dates)
  }
  node_ode_size <- length(ode_names)
  env_state <- patch$ode_state # fresh patch has no nodes: environment ODE only
  pr_surv0 <- patch$pr_survival(0)
  dens0 <- patch$density(0)

  n <- integer(n_spp)
  ode_chunks <- vector("list", n_spp)
  node_times <- vector("list", n_spp)
  patch_density <- vector("list", n_spp)
  pr_patch_survival <- vector("list", n_spp)
  for (i in seq_len(n_spp)) {
    h <- heights[[i]]
    ld <- log_densities[[i]]
    if (length(h) != length(ld)) {
      stop("heights and densities must have equal length (species ", i, ")")
    }
    ## The model requires nodes ordered by decreasing height.
    o <- order(h, decreasing = TRUE)
    h <- h[o]
    ld <- ld[o]
    mat <- matrix(0, nrow = node_ode_size, ncol = length(h))
    mat[hi, ] <- h
    times <- rep(0, length(h))
    if (birth_date) {
      times <- birth_dates[[i]][o]
      width <- diff(c(times, 0))
      if (length(times) != length(h) || any(width <= 0)) {
        stop("birth_dates must be before 0 and ascend as the heights descend ",
             "(species ", i, ")")
      }
      drivers <- patch$species[[i]]$extrinsic_drivers
      birth_rate <- vapply(times, function(t) drivers$evaluate("birth_rate", t),
                           numeric(1))
      mat[ldi, ] <- log(birth_rate) - ld
      mat[ii, ] <- rbind(width, width^2 / 2)
    } else {
      mat[ldi, ] <- ld
    }
    ode_chunks[[i]] <- as.vector(mat) # column-major: node-by-node, matching set_ode_state
    n[i] <- length(h)
    node_times[[i]] <- times
    patch_density[[i]] <- rep(dens0, length(h))
    pr_patch_survival[[i]] <- rep(pr_surv0, length(h))
  }

  list(
    time = 0,
    n = n,
    ode_state = c(unlist(ode_chunks, use.names = FALSE), env_state),
    node_times = node_times,
    patch_density = patch_density,
    pr_patch_survival = pr_patch_survival,
    ## Continue recruitment for t > 0; the seeded distribution replaces the t=0 recruit.
    node_schedule_times = lapply(p$node_schedule_times, function(tt) tt[tt > 1e-8])
  )
}

##' Write an exported patch state (from \code{\link{export_patch_state}}) into a
##' \code{Parameters} object so that the next \code{\link{run_scm}} starts from
##' that state instead of an empty patch. The state is carried on the
##' \code{Parameters} object (rather than passed separately) so the run stays
##' self-describing and reproducible, and so the seeding survives the reset at
##' the start of every run / schedule refinement.
##'
##' Each seeded cohort's patch density and survival at birth are read from the
##' resuming parameters' own regime at its birth date, not taken from the
##' export: a cohort seeded from a burn-in into a \code{\link{fitness_window}}
##' run weighs what the window gives it. Under the regime that exported it they
##' are the same numbers. The exported \code{patch_density} and
##' \code{pr_patch_survival} are checked only for their lengths.
##'
##' @title Seed Parameters with an initial patch state
##' @param p A \code{Parameters} object (its strategies must match the exported
##'   state).
##' @param state An exported state list from \code{\link{export_patch_state}}.
##' @return The modified \code{Parameters} object. The state's \code{events}
##'   are not carried on it: pass them to \code{\link{run_scm}} as
##'   \code{events = state$events}.
##' @seealso \code{\link{export_patch_state}}, \code{\link{run_scm}}
##' @export
set_initial_state <- function(p, state) {
  n_spp <- length(p$strategies)
  if (length(state$n) != n_spp) {
    stop("State has ", length(state$n),
         " species but Parameters has ", n_spp, " strategies")
  }
  p$initial_state <- state$ode_state
  p$n_initial_cohorts <- as.integer(state$n)
  p$initial_node_times <- unlist(state$node_times, use.names = FALSE)
  p$initial_patch_density <- unlist(state$patch_density, use.names = FALSE)
  p$initial_pr_patch_survival <- unlist(state$pr_patch_survival, use.names = FALSE)
  p$initial_time <- state$time
  if (!is.null(state$node_schedule_times)) {
    p$node_schedule_times <- state$node_schedule_times
  }
  p
}


##' Weight arrivals by a window of the run: the continuous patch's regime
##' (plant issue #519). Where the default \code{"meta-population"} patch type
##' weights a cohort by the density of patches at its age under a Weibull
##' disturbance regime, this weights it by when it arrived in one long run: a
##' cohort born inside \code{[start, end]} carries the tabulated weight at its
##' birth date and a cohort born outside carries zero. The net reproduction
##' ratio is then the weighted mean, over the cohorts the run introduced, of a
##' cohort's lifetime output: output times weight integrated over the birth
##' dates between the window's first and last cohorts (by establishment weight
##' on the birth-date coordinate), divided by the trapezium of the weight on
##' the same grid, so a uniform window is the exact mean however its edges sit
##' against the introduction grid. Patch survival is one.
##'
##' With no knots the window is uniform, each birth date weighted
##' \code{1 / (end - start)}. Knots give a piecewise-linear weight, held at the
##' end knots' values beyond them: a germination opportunity built from the
##' same rainfall series that drives the run, say.
##'
##' The run must extend past \code{end} by the longest lifespan that matters,
##' or the last arrivals' output is truncated; and the window should open after
##' the stand's burn-in, which \code{\link{export_patch_state}} lets a run pay
##' once.
##'
##' The window is closed at both ends, because the trapezium over
##' \code{[start, end]} needs a cohort at each. Invaders introduced only inside
##' the window need one at \code{end} too:
##' \code{node_schedule_times_periodic(start, end, by, include_end = TRUE)}.
##'
##' @title Weight arrivals by a window
##' @param p A \code{Parameters} object.
##' @param start,end The window, in years of the run.
##' @param weight_x,weight_y Optional knots of the weight over the window;
##'   \code{weight_x} increasing, \code{weight_y} non-negative. Omit both for a
##'   uniform window.
##' @return The modified \code{Parameters} object, with
##'   \code{patch_type = "fitness_window"}.
##' @seealso \code{\link{node_schedule_times_periodic}}
##' @export
fitness_window <- function(p, start, end, weight_x = numeric(0),
                          weight_y = numeric(0)) {
  if (!(is.numeric(start) && is.numeric(end) && length(start) == 1L &&
        length(end) == 1L && end > start)) {
    stop("`end` must exceed `start`")
  }
  if (length(weight_x) != length(weight_y)) {
    stop("`weight_x` and `weight_y` must have the same length")
  }
  p$patch_type <- "fitness_window"
  p$fitness_window_start <- as.numeric(start)
  p$fitness_window_end <- as.numeric(end)
  p$fitness_window_weight_x <- as.numeric(weight_x)
  p$fitness_window_weight_y <- as.numeric(weight_y)
  p
}

##' Node introduction times on a uniform grid, for a run that introduces a
##' cohort periodically rather than on the front-loaded default schedule of
##' \code{\link{node_schedule_times_default}}: the continuous patch (plant issue
##' #519). Times run from \code{start} in steps of \code{by}. By default they
##' stop short of \code{end}, so consecutive spans tile without repeating a
##' time: \code{periodic(0, 100, 1)} and \code{periodic(100, 200, 1)} meet at
##' 100 once. A cohort introduced at a run's \code{max_patch_lifetime} would
##' never be integrated either.
##'
##' A schedule for a \code{\link{fitness_window}} over \code{[start, end]} is
##' the other case: the window weighs a cohort at each end, so invaders
##' introduced only inside it need \code{include_end = TRUE} to have the
##' cohort at \code{end} the resident has. Without it their window average
##' misses the last arrival (2 \% low for a copy of the resident over ten).
##'
##' @title Periodic node introduction times
##' @param start First introduction time.
##' @param end The span's end.
##' @param by Spacing between introductions, in years.
##' @param include_end Whether a time landing on \code{end} is kept.
##' @return A numeric vector of introduction times, for one species of
##'   \code{p$node_schedule_times}.
##' @seealso \code{\link{fitness_window}}
##' @export
node_schedule_times_periodic <- function(start, end, by, include_end = FALSE) {
  if (!(by > 0)) {
    stop("`by` must be positive")
  }
  if (!(end > start)) {
    stop("`end` must exceed `start`")
  }
  times <- seq(start, end, by = by)
  if (include_end) times[times <= end] else times[times < end]
}
