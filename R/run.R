## One verb for running plant (#428). A run is described by a `scenario` -- the
## Parameters, Environment, Control and Events that go into it -- and returns a
## `plant_run`: what was run, what it reported, and the live solver for the
## operations that need it (invaders, gradients, export, continuing).

##' Describe a run of plant
##'
##' A scenario holds everything a run takes: who is present and how fitness is
##' weighted (\code{p}), the forcing (\code{env} and \code{events}) and the
##' numerics (\code{control}). \code{\link{run}} takes one, or builds one from
##' the same arguments, and every result carries the scenario it ran as
##' \code{$scenario}, so \code{run(res$scenario)} runs it again in the same
##' session. A scenario's environment does not survive \code{saveRDS()}.
##'
##' @param p A \code{Parameters} object.
##' @param env An \code{Environment}; by default the model's own.
##' @param control A \code{Control} object.
##' @param events An \code{\link{events}} object, or \code{NULL} for the
##'   introductions in \code{p$node_schedule_times} alone.
##' @return A list of class \code{plant_scenario} with \code{p}, \code{env},
##'   \code{control} and \code{events}.
##' @seealso \code{\link{run}}
##' @export
scenario <- function(p, env = NULL, control = plant::control(), events = NULL) {
  if (!inherits(p, "Parameters")) {
    stop("`p` must be a Parameters object", call. = FALSE)
  }
  model <- extract_RcppR6_template_types(p, "Parameters")[[1]]
  if (is.null(env)) {
    env <- Environment(model)
  }
  if (!inherits(control, "Control")) {
    stop("`control` must be a Control object", call. = FALSE)
  }
  if (!is.null(events) && !inherits(events, "Events")) {
    stop("`events` must be an Events object, as events() returns, or NULL",
         call. = FALSE)
  }
  structure(list(p = p, env = env, control = control, events = events),
            class = "plant_scenario")
}

##' @export
print.plant_scenario <- function(x, ...) {
  model <- extract_RcppR6_template_types(x$p, "Parameters")[[1]]
  n_events <- if (is.null(x$events)) 0L else length(x$events$time)
  cat(sprintf("<plant_scenario> %s, %d species, %g years, %d event(s)\n",
              model, length(x$p$strategies), x$p$max_patch_lifetime, n_events))
  invisible(x)
}

##' Run plant
##'
##' One verb for a run, whatever describes it:
##'
##' \itemize{
##'   \item a \code{\link{scenario}}, the bundle of everything a run takes;
##'   \item a \code{Parameters} object, with the rest given as arguments
##'     (\code{run(p, env = env)} is \code{run(scenario(p, env))});
##'   \item a previous result, which continues it to \code{until};
##'   \item an \code{Individual}, a single plant grown on its own.
##' }
##'
##' \code{method} chooses the solver: \code{"deterministic"}, by the method of
##' characteristics, or \code{"stochastic"}, the finite population of
##' individuals, arriving at random and each establishing and dying by chance.
##' The stochastic run takes the scenario's \code{p}, \code{env} and
##' \code{control} and runs the whole lifetime. Arrivals are drawn from each
##' species' birth rate and \code{p$patch_area}, not from
##' \code{p$node_schedule_times}; the result's \code{scenario$p} carries the
##' times drawn. It has no fitness integral, so
##' its result's fitness values are \code{NULL}, and it refuses \code{until},
##' \code{from}, \code{refine_schedule}, \code{record = "states"} and events.
##' Seed R's random number generator to repeat one.
##'
##' A single plant is grown by \code{run(individual, times = )}, which reports
##' it at each time, or \code{run(individual, sizes = )}, which grows it to
##' each size and bisects for the time it is reached. A plant grown on its own
##' sees the environment it is given, not a patch's shading.
##'
##' @param x What to run: a scenario, a \code{Parameters} object, a
##'   \code{plant_run} to continue, or an \code{Individual}.
##' @param ... For a \code{Parameters} object, the arguments below. For an
##'   \code{Individual} grown to \code{sizes}: \code{time_max}, the time to
##'   stop at if a size is never reached (default \code{Inf}); \code{warn},
##'   whether to warn when sizes are dropped (default \code{TRUE}); and
##'   \code{filter}, whether to drop sizes not reached from the result
##'   (default \code{FALSE}).
##' @param until Where to stop: \code{NULL} runs the whole patch lifetime, and a
##'   time stops at the first step boundary at or past it. A time at or before
##'   the start builds the run without taking a step, for a caller that drives
##'   it: \code{\link{run_until_converged}}, or \code{run(res, until = )}.
##'   Continue with \code{run(res)}.
##' @param from Where to start: \code{NULL} starts the patch empty; a state
##'   exported with \code{\link{export_patch_state}}, or a \code{plant_run}
##'   (exported where it stands), resumes from it. The resumed run introduces
##'   nodes on the scenario's own schedule from the export on, and applies the
##'   events the state still had to come, so no pulse or harvest pending at the
##'   export is lost. Those events are the run's: a scenario given with
##'   \code{from} must carry none of its own. The schedule is
##'   \code{p$node_schedule_times}, so a run whose introductions came from its
##'   own \code{events} instead resumes on \code{p}'s; set
##'   \code{p$node_schedule_times} to its times.
##' @param record What to keep beyond the outputs every run reports:
##'   \code{"steps"}, the patch at every step, which \code{tidy()} reads; and
##'   \code{"states"}, the ODE state at every step, which a gradient sweeps
##'   (\code{\link{stand_gradient}}).
##' @param method \code{"deterministic"} (the default) or \code{"stochastic"}.
##' @param env,control,events For a \code{Parameters} object, as
##'   \code{\link{scenario}} takes them. For an \code{Individual}, \code{env} is
##'   the \code{Environment} it grows in, by default its model's.
##' @param times,sizes For an \code{Individual}, the times to report it at, or
##'   the sizes to grow it to (increasing); give one.
##' @param size_name The size \code{sizes} are of.
##' @param refine_schedule Refine the node-introduction schedule before the
##'   run, bisecting the intervals whose nodes' error exceeds
##'   \code{control$schedule_eps}, up to \code{control$schedule_nsteps} times.
##'   The refined schedule is in the result's \code{scenario$p}, with the ODE
##'   schedule the final run took (\code{ode_times}, \code{ode_step_sizes},
##'   \code{ode_soil_input_slopes}, \code{ode_soil_substep_ends}), so a later
##'   run of those parameters takes those steps exactly. With
##'   \code{record = "states"}, refinement's own runs keep no states and one
##'   more run keeps them.
##' @return A \code{plant_run}: a list with \code{scenario} (what was run, its
##'   \code{p} the solver's, carrying any refined schedule), \code{time} and
##'   \code{complete} (where the run stands), \code{offspring_production} and
##'   \code{net_reproduction_ratios} (per species, at \code{time}),
##'   \code{event_log} (what each event did), \code{mutants} (the invaders'
##'   \code{Parameters} for a \code{\link{run_mutant}} result, else \code{NULL}),
##'   \code{method}, and \code{solver}, the live \code{SCM}. A stochastic
##'   result's fitness values and \code{event_log} are \code{NULL}, its
##'   \code{steps} hold what \code{record = "steps"} kept, and its solver is a
##'   \code{StochasticPatchRunner}. An \code{Individual}'s run returns a
##'   \code{plant_individual_run} instead: a list with \code{time} (the times given, or
##'   when each size was reached, \code{NA} for one that was not), \code{state} (the ODE
##'   state at each time, as a matrix), \code{individual} (a list of the plant
##'   at each) and \code{env}, and when grown to sizes \code{trajectory}, the
##'   state at every step taken.
##'
##'   A result holds until its solver moves. Whatever moves it -- \code{run()}
##'   continuing it, \code{run_mutant()}, \code{run_until_converged()} --
##'   returns the result that describes it, and every function that takes a
##'   result refuses one whose solver has moved on since.
##'   \code{fixed_point_fitness()} uses a result up, so give it every invader
##'   in one trait matrix. The values other than
##'   \code{solver} are copies, so \code{saveRDS()} keeps them; the scenario's
##'   environment and the solver cannot be used after reading back.
##' @seealso \code{\link{scenario}}, \code{\link{run_mutant}},
##'   \code{\link{export_patch_state}}, \code{\link{stand_gradient}}
##' @export
run <- function(x, ...) {
  UseMethod("run")
}

##' @rdname run
##' @export
run.plant_scenario <- function(x, ..., until = NULL, from = NULL,
                               record = character(0), refine_schedule = FALSE,
                               method = c("deterministic", "stochastic")) {
  no_more_arguments(...)
  method <- match.arg(method)
  unknown <- setdiff(record, c("steps", "states"))
  if (length(unknown) > 0L) {
    stop("`record` takes \"steps\" and \"states\", not ",
         paste0("\"", unknown, "\"", collapse = ", "), call. = FALSE)
  }
  if (method == "stochastic") {
    return(run_stochastic(x, until, from, record, refine_schedule))
  }
  if (refine_schedule && !is.null(until)) {
    stop("`refine_schedule` refines over the whole run, so it cannot stop at `until`",
         call. = FALSE)
  }
  check_until(until)
  if (!is.null(from)) {
    x <- resume_scenario(x, from)
  }

  types <- extract_RcppR6_template_types(x$p, "Parameters")
  events <- if (is.null(x$events)) empty_events() else x$events
  scm <- do.call("SCM", types)(x$p, x$env, events, x$control)
  scm$collect <- "steps" %in% record
  keep_states <- "states" %in% record
  if (refine_schedule) {
    ## Refinement's runs are bisected against and discarded, so they keep no
    ## states; the run a sweep walks is one more, after the schedule settles.
    scm$refine_schedule()
    if (keep_states) {
      scm$keep_states <- TRUE
      scm$run()
    }
  } else {
    scm$keep_states <- keep_states
    if (is.null(until)) {
      scm$run()
    } else if (until > x$p$initial_time) {
      scm$run_to(until)
    }
  }
  x$p <- scm$parameters
  new_plant_run(scm, x, record)
}

##' @rdname run
##' @export
run.Parameters <- function(x, env = NULL, control = plant::control(), events = NULL,
                           ...) {
  run(scenario(x, env = env, control = control, events = events), ...)
}

##' @rdname run
##' @export
run.plant_run <- function(x, ..., until = NULL) {
  if (...length() > 0L) {
    stop("A run is continued with `until` alone; to change anything else, ",
         "run its scenario again", call. = FALSE)
  }
  check_until(until)
  scm <- as_scm(x)
  if (!is.null(until) && until < x$time) {
    stop("This run stands at t = ", format(x$time), ", past `until`", call. = FALSE)
  }
  if (is.null(until) || until > x$time) {
    scm$run_to(if (is.null(until)) Inf else until)
  }
  new_plant_run(scm, x$scenario, x$record, x$mutants)
}

## The scenario a resume runs: the state on `p`, introductions on the
## scenario's own schedule from the export on, and the state's pending events.
resume_scenario <- function(x, from) {
  if (!is.null(x$events)) {
    stop("With `from`, the run's events are the state's pending ones; give the ",
         "scenario none of its own", call. = FALSE)
  }
  if (inherits(from, "plant_run")) {
    from <- export_patch_state(from)
  }
  if (!is.list(from) || is.null(from$ode_state) || inherits(from, "R6")) {
    stop("`from` must be a state from export_patch_state() or make_initial_state(), ",
         "or a plant_run", call. = FALSE)
  }
  schedule <- x$p$node_schedule_times
  x$p <- set_initial_state(x$p, from)
  ## The export's residual schedule is the exporting run's; the resumed run
  ## introduces on its own, so a scenario longer than the exporting run keeps
  ## recruiting. On the exporting scenario the two are the same times.
  if (length(schedule) > 0L) {
    x$p$node_schedule_times <- Map(function(tt, made) {
      tt[tt >= from$time - 1e-8 & !(tt %in% made)]
    }, schedule, from$node_times)
  }
  pending <- from$events
  if (length(pending$time) > 0L) {
    ## An events list replaces p's schedule as the source of introductions, so
    ## the introductions in it are rebuilt from the schedule just set.
    keep <- pending$type != "node_introduction"
    actions <- Events(time = pending$time[keep], type = pending$type[keep],
                      target = pending$target[keep],
                      target_index = pending$target_index[keep],
                      params = pending$params[keep])
    x$events <- events(node_introductions(x$p), actions)
  }
  x
}

check_until <- function(until) {
  if (!is.null(until) && !(is.numeric(until) && length(until) == 1L && !is.na(until))) {
    stop("`until` must be one time, or NULL for the whole run", call. = FALSE)
  }
}

no_more_arguments <- function(...) {
  if (...length() > 0L) {
    nms <- names(list(...))
    if (is.null(nms)) nms <- rep("", ...length())
    nms[nms == ""] <- "an unnamed argument"
    stop("run() does not take ", paste(nms, collapse = ", "), " here", call. = FALSE)
  }
}

## The stochastic solver under run(). It runs the whole lifetime and reports no
## fitness, so whatever needs a stopping point, a resume, a recording for
## gradients or events is refused by name.
run_stochastic <- function(x, until, from, record, refine_schedule) {
  refused <- c(until = !is.null(until), from = !is.null(from),
               refine_schedule = isTRUE(refine_schedule),
               `record = "states"` = "states" %in% record,
               events = !is.null(x$events))
  if (any(refused)) {
    stop("The stochastic solver runs the whole lifetime with no events; it ",
         "cannot take ", paste(names(refused)[refused], collapse = ", "),
         call. = FALSE)
  }
  obj <- stochastic_runner(x$p, x$env, x$control)
  steps <- NULL
  if ("steps" %in% record) {
    steps <- collect_stochastic_steps(obj)
  } else {
    obj$run()
  }
  ## The arrivals drawn are on the solver's parameters; keep them, so a saved
  ## result records who arrived when.
  x$p <- obj$parameters
  structure(list(scenario = x, time = obj$time, complete = obj$complete,
                 offspring_production = NULL, net_reproduction_ratios = NULL,
                 event_log = NULL, record = record, mutants = NULL,
                 method = "stochastic", steps = steps, solver = obj),
            class = "plant_run")
}

new_plant_run <- function(scm, scenario, record, mutants = NULL) {
  patch <- scm$patch
  structure(list(scenario = scenario,
                 time = scm$time,
                 complete = scm$complete,
                 offspring_production = patch$offspring_production,
                 net_reproduction_ratios = patch$net_reproduction_ratios,
                 event_log = scm$event_log,
                 record = record,
                 mutants = mutants,
                 method = "deterministic",
                 solver = scm),
            class = "plant_run")
}

##' @rdname run
##' @export
run.Individual <- function(x, env = NULL, ..., times = NULL, sizes = NULL,
                           size_name = "height") {
  if (is.null(times) == is.null(sizes)) {
    stop("Give one of `times` and `sizes`", call. = FALSE)
  }
  if (is.null(env)) {
    env <- Environment(x$strategy_name)
  }
  out <- if (is.null(times)) {
    individual_grow_to_sizes(x, sizes, size_name, env, ...)
  } else {
    no_more_arguments(...)
    individual_grow_to_times(x, times, env)
  }
  structure(out, class = "plant_individual_run")
}

##' @export
print.plant_run <- function(x, ...) {
  model <- extract_RcppR6_template_types(x$scenario$p, "Parameters")[[1]]
  cat(sprintf("<plant_run> %s%s, %d species%s, at t = %g of %g%s\n", model,
              if (identical(x$method, "stochastic")) " (stochastic)" else "",
              length((x$mutants %||% x$scenario$p)$strategies),
              if (is.null(x$mutants)) "" else " invading",
              x$time, x$scenario$p$max_patch_lifetime,
              if (isTRUE(x$complete)) " (complete)" else ""))
  if (!is.null(x$offspring_production)) {
    cat("  offspring_production:", format(x$offspring_production, digits = 4), "\n")
  }
  invisible(x)
}

##' @importFrom generics tidy
##' @export
generics::tidy

##' Tidy a run's recorded steps
##'
##' @param x A \code{plant_run} run with \code{record = "steps"}, or a single
##'   plant's run.
##' @param ... Unused.
##' @return For a deterministic run, a list: the tidied patch at every step,
##'   \code{steps}, \code{n_spp}, \code{species} and \code{env} (as
##'   \code{\link{tidy_patch}} gives them); then the run's
##'   \code{offspring_production} and \code{net_reproduction_ratios}; \code{p},
##'   the solver's parameters (carrying any refined schedule); \code{events},
##'   the events it was given; and \code{event_log}, what each event did.
##'   For a stochastic run, \code{time}, \code{species} (an array per species,
##'   \code{[variable, time, individual]}, with an \code{is_alive} attribute)
##'   and \code{env}. For a single plant, a tibble of its state at each
##'   reported time, as \code{\link{tidy_individual}} gives it.
##' @export
tidy.plant_run <- function(x, ...) {
  ## A deterministic result is checked as current first: its steps are the
  ## solver's history, which a later run overwrites.
  if (!identical(x$method, "stochastic")) {
    as_scm(x)
  }
  if (!"steps" %in% x$record) {
    stop("This run kept no steps: run it with record = \"steps\"", call. = FALSE)
  }
  if (identical(x$method, "stochastic")) {
    return(tidy_stochastic_steps(x$steps))
  }
  scm <- x$solver
  if (length(scm$history) == 0L) {
    stop("This run has taken no step yet: continue it with run(res) first",
         call. = FALSE)
  }
  results <- tidy_patch(lapply(scm$history, "[[", "state"))
  patch <- scm$patch
  results[["offspring_production"]] <- patch$offspring_production
  results[["net_reproduction_ratios"]] <- patch$net_reproduction_ratios
  results[["p"]] <- scm$parameters
  ## Events are given separately from `p`, so `p` alone does not describe the
  ## run: the schedule asked for and what was applied (including how much of
  ## each pulse the soil accepted and shed) travel with the tables.
  results[["events"]] <- scm$events
  results[["event_log"]] <- scm$event_log
  results
}

##' @rdname tidy.plant_run
##' @export
tidy.plant_individual_run <- function(x, ...) {
  tidy_individual(x)
}

##' Run invaders against a resident
##'
##' The invaders in \code{p} run in the field the resident recorded, which they
##' read and do not move; with \code{constant = TRUE}, in the one field the
##' resident stands in now (see \code{\link{fixed_point_fitness}}).
##'
##' @param resident A complete \code{plant_run} of the resident.
##' @param p The invaders' \code{Parameters}, as \code{\link{add_mutant}} builds
##'   them.
##' @param constant Hold the field the resident stands in now, rather than
##'   replay the one it recorded.
##' @return A \code{plant_run} of the invaders: their outputs, with
##'   \code{mutants} set to \code{p} and \code{scenario} the resident's. It
##'   shares the resident's solver, which now holds the invaders, so the
##'   resident's result is refused from here on (see \code{\link{run}}).
##' @export
run_mutant <- function(resident, p, constant = FALSE) {
  if (!inherits(resident, "plant_run") || !is.null(resident$mutants)) {
    stop("`resident` must be a resident's plant_run, as run() returns", call. = FALSE)
  }
  scm <- as_scm(resident)
  if (!isTRUE(resident$complete)) {
    stop("The resident has not finished its run; invaders read its whole recording",
         call. = FALSE)
  }
  if (constant) {
    scm$run_mutant_constant(p)
  } else {
    scm$run_mutant(p)
  }
  new_plant_run(scm, resident$scenario, resident$record, mutants = p)
}

## The solver behind a run, for the functions that take either. A result is
## refused once its solver has moved on, since everything read from it would
## describe some later run.
as_scm <- function(x) {
  if (!inherits(x, "plant_run")) {
    return(x)
  }
  if (identical(x$method, "stochastic")) {
    stop("This is a stochastic run, which has no fitness, gradient, export or continuation: ",
         "that needs the deterministic solver", call. = FALSE)
  }
  scm <- x$solver
  time <- tryCatch(scm$time, error = function(e) {
    stop("This result's solver is gone (a result read back with readRDS() keeps ",
         "its values only): run its scenario again", call. = FALSE)
  })
  if (!identical(time, x$time) ||
      !identical(scm$parameters, x$mutants %||% x$scenario$p)) {
    stop("This result's solver has moved on since it was returned (continued, ",
         "converged or invaded): use the result of whatever moved it", call. = FALSE)
  }
  scm
}
