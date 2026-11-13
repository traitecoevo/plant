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
##'   \item a previous result, which continues it to \code{until}.
##' }
##'
##' The run is the deterministic solver, by the method of characteristics.
##'
##' @param x What to run: a scenario, a \code{Parameters} object, or a
##'   \code{plant_run} to continue.
##' @param ... For a \code{Parameters} object, \code{env}, \code{control} and
##'   \code{events}, as \code{\link{scenario}} takes them, then the arguments
##'   below.
##' @param until Where to stop: \code{NULL} runs the whole patch lifetime, and a
##'   time stops at the first step boundary at or past it. A time at or before
##'   the start builds the run without taking a step, as
##'   \code{run_scm(run = FALSE)} does. Continue with \code{run(res)}.
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
##' @param refine_schedule Refine the node-introduction schedule before the
##'   run, as \code{\link{run_scm}} does. The refined schedule is in the result's
##'   \code{scenario$p}.
##' @return A \code{plant_run}: a list with \code{scenario} (what was run, its
##'   \code{p} the solver's, carrying any refined schedule), \code{time} and
##'   \code{complete} (where the run stands), \code{offspring_production} and
##'   \code{net_reproduction_ratios} (per species, at \code{time}),
##'   \code{event_log} (what each event did), \code{mutants} (the invaders'
##'   \code{Parameters} for a \code{\link{run_mutant}} result, else \code{NULL})
##'   and \code{solver}, the live \code{SCM}.
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
                               record = character(0), refine_schedule = FALSE) {
  no_more_arguments(...)
  unknown <- setdiff(record, c("steps", "states"))
  if (length(unknown) > 0L) {
    stop("`record` takes \"steps\" and \"states\", not ",
         paste0("\"", unknown, "\"", collapse = ", "), call. = FALSE)
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
    ## As run_scm(): refinement's runs keep no states; the run a sweep walks is
    ## one more, after the schedule settles.
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

new_plant_run <- function(scm, scenario, record, mutants = NULL) {
  structure(list(scenario = scenario,
                 time = scm$time,
                 complete = scm$complete,
                 offspring_production = scm$offspring_production,
                 net_reproduction_ratios = scm$net_reproduction_ratios,
                 event_log = scm$event_log,
                 record = record,
                 mutants = mutants,
                 solver = scm),
            class = "plant_run")
}

##' @export
print.plant_run <- function(x, ...) {
  model <- extract_RcppR6_template_types(x$scenario$p, "Parameters")[[1]]
  cat(sprintf("<plant_run> %s, %d species%s, at t = %g of %g%s\n", model,
              length(x$offspring_production),
              if (is.null(x$mutants)) "" else " invading",
              x$time, x$scenario$p$max_patch_lifetime,
              if (isTRUE(x$complete)) " (complete)" else ""))
  cat("  offspring_production:", format(x$offspring_production, digits = 4), "\n")
  invisible(x)
}

##' @importFrom generics tidy
##' @export
generics::tidy

##' Tidy a run's recorded steps
##'
##' @param x A \code{plant_run} run with \code{record = "steps"}.
##' @param ... Unused.
##' @return The tidied patch: \code{steps}, \code{n_spp}, \code{species} and
##'   \code{env}, as \code{run_scm(collect = TRUE)} gives them.
##' @export
tidy.plant_run <- function(x, ...) {
  scm <- as_scm(x)
  if (!"steps" %in% x$record) {
    stop("This run kept no steps: run it with record = \"steps\"", call. = FALSE)
  }
  tidy_patch(lapply(scm$history, "[[", "state"))
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
