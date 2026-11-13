##' Construct a \code{Control} object. \code{control()} is a lowercase alias for
##' the \code{Control()} constructor, whose defaults are the pragmatic,
##' fast-ish settings used for essentially all of plant's runs (see
##' \code{control.cpp}). \code{control_accurate()} tightens the ODE and schedule
##' tolerances for high-accuracy runs at the cost of speed. \code{control_tf24()}
##' sets TF24's step control: the relative and absolute tolerances at
##' \code{tol}, the soil substepped, and the soil layers' tolerance at 10 times
##' that. Against \code{Control()} on a TF24 stand under daily rain it takes
##' about 60\% of the steps, at about a third of the error.
##'
##' \code{ode_tol_scale_soil} multiplies the tolerance on TF24's soil layers in
##' the SCM's adaptive steps, which decides which steps it takes and nothing
##' else. Its default of 1 leaves every step as it was. The stochastic runner
##' ignores it, and warns. \code{ode_soil_substep_tol} is the relative
##' tolerance that chooses the soil's substeps inside a step (default 1e-9).
##'
##' @title Control presets
##' @param ... Named control fields, passed to \code{Control()}.
##' @param base An optional \code{Control} object to start from; defaults are
##'   used if omitted.
##' @return A \code{Control} object.
##' @rdname control_presets
##' @export
control <- function(...) Control(...)

##' @rdname control_presets
##' @export
control_accurate <- function(base = Control()) {
  base$ode_tol_rel       <- 1e-7
  base$ode_tol_abs       <- 1e-7
  base$ode_step_size_max <- 1e-1
  base$schedule_eps      <- 1e-3
  base
}

##' @rdname control_presets
##' @param tol The relative tolerance, \code{ode_tol_rel}.
##' @export
control_tf24 <- function(tol = 1e-5, base = Control()) {
  base$ode_tol_rel         <- tol
  base$ode_tol_abs         <- tol
  base$ode_soil_substep    <- TRUE
  base$ode_tol_scale_soil  <- 10
  base
}


##' A model's bare \code{Parameters}: its defaults, with no strategies and no
##' trade-offs attached, so a trait given to \code{\link{add_strategies}} is
##' written straight to \code{pars}. \code{\link{default_parameters}} attaches
##' plant's testing trade-offs to the same object. The environment named here
##' only picks the template; a run takes its \code{Environment} separately (see
##' \code{\link{run}}).
##' @title A model's bare parameters
##' @author Rich FitzJohn
##' @param type Any model name as a string, e.g.: \code{"FF16"}.
##' @param env The environment type's name; by default the model's own.
##' @return A \code{Parameters} object.
##' @export
base_parameters <- function(type = NA, env = environment_type(type)) {
  Parameters(type, env)()
}


##' Export the full state of a patch from a (run) \code{SCM} so it can be
##' re-imported to initialise a new run (see \code{\link{set_initial_state}}). The
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
##' @param scm A \code{plant_run}, or an \code{SCM} object that has been run.
##' @param step Optional 1-based index into \code{scm$history} (requires the run
##'   to have recorded its steps, \code{record = "steps"}); when \code{NULL}
##'   (default) the SCM's current/final patch is used.
##' @return A list describing the patch state, suitable for
##'   \code{\link{set_initial_state}}: \code{time}, \code{n} (nodes per species),
##'   \code{ode_state} (flat), the per-species lists \code{node_times},
##'   \code{patch_density} and \code{pr_patch_survival}, the residual
##'   \code{node_schedule_times}, and \code{events}, the run's events still to
##'   come -- its residual introductions and every action the run had not yet
##'   applied when the state was taken (an action at the export time itself is
##'   pending in a collected snapshot and in a run that \code{SCM$run_to()}
##'   stopped there, and applied in a finished run's final patch). Events are a \code{\link{run}} argument rather than a field
##'   of \code{Parameters}: \code{run(p, from = state)} applies them, and a
##'   run of \code{set_initial_state(p, state)} must be given
##'   \code{events = state$events}, or it re-runs the residual introductions
##'   alone and silently loses any pulse or harvest still ahead.
##'
##' With the absorbing top layer on (\code{Control$absorb_epsilon}), each
##' species' block of \code{ode_state} ends with the layer: its pinned
##' individual's states, then its density, its two fitness weights and its two
##' seed totals. A state exported without the layer has no such block, so
##' resuming it with the layer on needs that many zeros inserted after each
##' species' nodes.
##' @seealso \code{\link{set_initial_state}}, \code{\link{run}}
##' @export
export_patch_state <- function(scm, step = NULL) {
  scm <- as_scm(scm)
  patch <- if (is.null(step)) scm$patch else scm$history[[step]]
  time <- patch$time
  species <- patch$species
  tol <- 1e-8

  ## Residual schedule: the introductions not yet made. An initialised node carries
  ## its introduction time, so the residual is every scheduled time at or after
  ## `time` that no initialised node of that species carries. A patch snapshot sits
  ## integrated *to* its next introduction time without having introduced that
  ## node yet, so the node due at `time` itself stays in the residual; a run
  ## that ended on an introduction (one scheduled at max_patch_lifetime) has
  ## made it, and a time an initialised node carries is never residual.
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

##' Build an initial patch state from a specified size distribution, for initialising
##' a patch at age 0 with pre-existing plants instead of growing it from empty
##' (the ecological motivation of \code{plant} issue #304). The returned object
##' is consumed by \code{\link{set_initial_state}}.
##'
##' Each species' initial nodes are described by their heights and densities;
##' the remaining ODE state (mortality, fecundity, accumulated reproduction,
##' heartwood, ...) starts at zero, all nodes are introduced at patch age 0, and
##' the recruitment schedule continues for \code{t > 0} (the initialised distribution
##' replaces the \code{t = 0} recruit). Pathologically large/dense initial
##' conditions can produce non-finite densities; \code{\link{run}} guards against
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
##' @return A state list suitable for \code{\link{set_initial_state}}.
##' @seealso \code{\link{set_initial_state}}, \code{\link{export_patch_state}}
##' @export
make_initial_state <- function(p, heights, densities = NULL,
                               log_densities = NULL, env = NULL,
                               ctrl = control()) {
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
  ldi <- match("log_density", ode_names)
  if (is.na(hi) || is.na(ldi)) {
    stop("could not locate 'height'/'log_density' in node ODE names")
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
    mat[ldi, ] <- ld
    ode_chunks[[i]] <- as.vector(mat) # column-major: node-by-node, matching set_ode_state
    n[i] <- length(h)
    node_times[[i]] <- rep(0, length(h))
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
    ## Continue recruitment for t > 0; the initialised distribution replaces the t=0 recruit.
    node_schedule_times = lapply(p$node_schedule_times, function(tt) tt[tt > 1e-8])
  )
}

##' Write an exported patch state (from \code{\link{export_patch_state}}) into a
##' \code{Parameters} object so that the next \code{\link{run}} starts from
##' that state instead of an empty patch. The state is carried on the
##' \code{Parameters} object (rather than passed separately) so the run stays
##' self-describing and reproducible, and so the initial state survives the reset at
##' the start of every run / schedule refinement.
##'
##' Each initial cohort's patch density and survival at birth are read from the
##' resuming parameters' own regime at its birth date, not taken from the
##' export: a cohort initialised from a burn-in into a \code{\link{fitness_window}}
##' run weighs what the window gives it. Under the regime that exported it they
##' are the same numbers. The exported \code{patch_density} and
##' \code{pr_patch_survival} are checked only for their lengths.
##'
##' @title Give Parameters an initial patch state
##' @param p A \code{Parameters} object (its strategies must match the exported
##'   state).
##' @param state An exported state list from \code{\link{export_patch_state}}.
##' @return The modified \code{Parameters} object. The state's \code{events}
##'   are not carried on it: pass them to \code{\link{run}} as
##'   \code{events = state$events}, or resume with \code{run(p, from = state)},
##'   which does both.
##' @seealso \code{\link{export_patch_state}}, \code{\link{run}}
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
##' (plant issue #519). Where the default \code{fitness_weighting = "patch_age"}
##' weights a cohort by the density of patches at its age under a Weibull
##' disturbance regime, this weights it by when it arrived in one long run: a
##' cohort born inside \code{[start, end]} carries the tabulated weight at its
##' birth date and a cohort born outside carries zero. The net reproduction
##' ratio is then the weighted mean, over the cohorts the run introduced, of a
##' cohort's lifetime output: the trapezium over birth dates of output times
##' weight, divided by the trapezium of the weight on the same grid, so a
##' uniform window is the exact mean however its edges sit against the
##' introduction grid. Patch survival is one.
##'
##' With no knots the window is uniform, each birth date weighted
##' \code{1 / (end - start)}. Knots give a piecewise-linear weight, constant at the
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
##'   \code{fitness_weighting = "window"}.
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
  p$fitness_weighting <- "window"
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
