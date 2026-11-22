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
##' \code{$scenario}, so \code{run(res$scenario)} runs it again.
##'
##' A scenario is a value: it keeps its environment as a recipe, plain values
##' for everything a run reads from it, and each run builds the environment
##' afresh from that recipe. So a scenario, and a result carrying one, survives
##' \code{saveRDS()} and runs after \code{readRDS()}, and a change made to
##' \code{env} after the scenario was made does not reach it.
##' \code{Environment(sc$env)} builds the environment a scenario describes.
##'
##' @param p A \code{Parameters} object.
##' @param env An \code{Environment}, or a recipe of one (a scenario's
##'   \code{env}); by default the model's own.
##' @param control A \code{Control} object.
##' @param events An \code{\link{events}} object, or \code{NULL} for the
##'   introductions in \code{p$node_schedule_times} alone.
##' @return A list of class \code{plant_scenario} with \code{p}, \code{env}
##'   (the environment's recipe), \code{control} and \code{events}.
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
  env <- environment_recipe(env)
  if (!identical(env$type, environment_type(model))) {
    stop("`env` is a ", env$type, " environment, and ", model, " runs in a ",
         environment_type(model), " one", call. = FALSE)
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
##'   \item an \code{Individual}, a single plant grown on its own;
##'   \item a \code{\link{landscape}}, several patches and the flows that link
##'     them.
##' }
##'
##' \code{method} chooses the solver: \code{"deterministic"}, by the method of
##' characteristics, or \code{"stochastic"}, the finite population of
##' individuals, arriving at random and each establishing and dying by chance.
##' The stochastic run takes the scenario's \code{p}, \code{env} and
##' \code{control} and runs the whole lifetime. Arrivals are drawn from each
##' species' birth rate (the curve the deterministic run's inflow reads) and
##' \code{p$patch_area}, as a Poisson process in continuous time, not from
##' \code{p$node_schedule_times}; the result's \code{scenario$p} carries the
##' times drawn. Introductions named in \code{events} replace them, one arrival
##' per time. It has no fitness integral, so
##' its result's fitness values are \code{NULL}, and it refuses \code{until},
##' \code{from}, \code{refine_schedule} and \code{record = "states"}. It takes
##' \code{events}' harvests, resource pulses and climate extremes, drawing each
##' individual's removal by chance. Its deaths are drawn every
##' \code{control$stochastic_death_interval} years, 0.01 by default, which is
##' set in the run's \code{control}, as in
##' \code{control = control(stochastic_death_interval = 0.05)}. Survival is
##' exact at any interval, but a plant that has died shades the living until
##' the next draw. Where shade is what kills, a longer interval runs faster
##' and leaves fewer plants than continuous checks would: in a closing FF16
##' canopy about 2\% fewer at 0.05 and 0.5\% at 0.01.
##' Seed R's random number generator to repeat one.
##'
##' A single plant is grown by \code{run(individual, times = )}, which reports
##' it at each time, or \code{run(individual, sizes = )}, which grows it to
##' each size: the time a size is reached is located on the step that crosses
##' it, and the run lands there, so the time, state and plant reported for a
##' size describe the same instant. A plant grown on its own
##' sees the environment it is given, not a patch's shading.
##'
##' A landscape is run patch by patch. Each patch is run whole, as its own
##' scenario would be, after the patches that feed it; what those passed down
##' their flows is added to its events. \code{record},
##' \code{refine_schedule} and \code{method} apply to every patch, and
##' \code{until} and \code{from} are refused.
##'
##' @param x What to run: a scenario, a \code{Parameters} object, a
##'   \code{plant_run} to continue, an \code{Individual}, or a landscape.
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
##'   the \code{Environment} it grows in, or a recipe of one, by default its
##'   model's.
##' @param times,sizes For an \code{Individual}, the times to report it at, or
##'   the sizes to grow it to (increasing); give one.
##' @param size_name The size \code{sizes} are of.
##' @param refine_schedule Refine the node-introduction schedule before the
##'   run, bisecting the intervals whose nodes' error exceeds
##'   \code{control$schedule_eps}, up to \code{control$schedule_nsteps} times.
##'   The result is a run of the refined schedule: if the last pass still
##'   bisected, that schedule is run once more. The schedule is in the
##'   result's \code{scenario$p}, with the ODE schedule that run took
##'   (\code{ode_times}, \code{ode_step_sizes}, \code{ode_soil_input_slopes},
##'   \code{ode_soil_substep_ends}), so a later run of those parameters takes
##'   those steps exactly. With \code{record = "states"}, refinement's own
##'   passes keep no states, and the run the result reports keeps them; if the
##'   refinement converged, that is one more run.
##' @return A \code{plant_run}: a list with \code{scenario} (what was run, its
##'   \code{p} the solver's, carrying any refined schedule), \code{time} and
##'   \code{complete} (where the run stands), \code{offspring_production} and
##'   \code{net_reproduction_ratios} (per species, at \code{time}),
##'   \code{event_log} (what each event did), \code{mutants} (the invaders'
##'   \code{Parameters} for a \code{\link{run_mutant}} result, else \code{NULL}),
##'   \code{method}, and \code{solver}, the live
##'   \code{DeterministicPatchRunner}. A stochastic
##'   result's fitness values are \code{NULL}, its \code{event_log} holds what
##'   its events did (a removal reports individuals drawn for, and those removed
##'   per unit area), its
##'   \code{steps} hold what \code{record = "steps"} kept, and its solver is a
##'   \code{StochasticPatchRunner}. An \code{Individual}'s run returns a
##'   \code{plant_individual_run} instead: a list with \code{time} (the times given, or
##'   when each size was reached, \code{NA} for one that was not), \code{state} (the ODE
##'   state at each time, as a matrix), \code{individual} (a list of the plant
##'   at each) and \code{env}, and when grown to sizes \code{trajectory}, the
##'   time and state at every step taken, the steps landed on each size
##'   among them.
##'
##'   A landscape's run returns a \code{plant_landscape_run}: a list with
##'   \code{landscape} (as given, so \code{run(res$landscape)} runs it again),
##'   \code{patches} (a \code{plant_run} per patch, by name; a patch's
##'   \code{scenario} is what that patch ran, its events including what
##'   reached it down the flows), \code{flows} (a data frame of every amount
##'   passed between patches: \code{time}, \code{from}, \code{to},
##'   \code{type}, \code{out_of}, \code{into} and \code{amount}, the amount
##'   that reached \code{to}) and \code{method}.
##'   \code{tidy()} of it is a list of each patch's.
##'
##'   A result holds until its solver moves. Whatever moves it -- \code{run()}
##'   continuing it, \code{run_mutant()}, \code{run_until_converged()} --
##'   returns the result that describes it, and every function that takes a
##'   result refuses one whose solver has moved on since.
##'   \code{fixed_point_fitness()} uses a result up, so give it every invader
##'   in one trait matrix.
##'
##'   \code{saveRDS()} keeps everything but the solver: the values are copies
##'   and the scenario holds its environment as a recipe. A deterministic
##'   result read back runs again from its scenario: \code{run(res, until =
##'   res$time)} rebuilds it where it stood, taking the same steps, and
##'   \code{run(res)} carries on to the end. That repeats the run's work. An
##'   invaders' result is rebuilt by running its resident again and then
##'   \code{run_mutant()}.
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
  solver <- do.call("DeterministicPatchRunner", types)(
    x$p, live_environment(x$env), events, x$control)
  solver$collect <- "steps" %in% record
  solver$keep_states <- "states" %in% record
  if (refine_schedule) {
    solver$refine_schedule()
  } else if (is.null(until)) {
    solver$run()
  } else if (until > x$p$initial_time) {
    solver$run_to(until)
  }
  x$p <- solver$parameters
  new_plant_run(solver, x, record)
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
  if (solver_is_gone(x)) {
    return(run_again(x, until))
  }
  solver <- as_solver(x)
  check_until_not_past(x, until)
  if (is.null(until) || until > x$time) {
    solver$run_to(if (is.null(until)) Inf else until)
  }
  new_plant_run(solver, x$scenario, x$record, x$mutants)
}

check_until_not_past <- function(x, until) {
  if (!is.null(until) && until < x$time) {
    stop("This run stands at t = ", format(x$time), ", past `until`", call. = FALSE)
  }
}

## A deterministic result whose solver is gone: one read back with readRDS().
solver_is_gone <- function(x) {
  !identical(x$method, "stochastic") &&
    !tryCatch({x$solver$time; TRUE}, error = function(e) FALSE)
}

## A result read back, run again from its scenario. Its scenario is the one the
## solver ran -- its `p` carries any refined schedule and the ODE steps that
## schedule took, and a resumed run's state and pending events -- so running it
## to the same stopping point takes the same steps and lands on the same step
## boundary, and running on past it is the run continued.
run_again <- function(x, until) {
  if (!is.null(x$mutants)) {
    stop("This invaders' result was read back without its solver: run its ",
         "resident's scenario again, then run_mutant()", call. = FALSE)
  }
  check_until_not_past(x, until)
  run(x$scenario, until = until, record = x$record)
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
    ## An events list naming introductions replaces p's schedule as their
    ## source, so the introductions in it are rebuilt from the schedule just set.
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
## fitness, so whatever needs a stopping point, a resume or a recording for
## gradients is refused by name. It takes a run's events as the deterministic
## solver does: introductions named there replace the drawn arrivals.
run_stochastic <- function(x, until, from, record, refine_schedule) {
  refused <- c(until = !is.null(until), from = !is.null(from),
               refine_schedule = isTRUE(refine_schedule),
               `record = "states"` = "states" %in% record)
  if (any(refused)) {
    stop("The stochastic solver runs the whole lifetime; it cannot take ",
         paste(names(refused)[refused], collapse = ", "), call. = FALSE)
  }
  obj <- stochastic_runner(x$p, live_environment(x$env), x$control, x$events)
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
                 event_log = obj$event_log, record = record, mutants = NULL,
                 method = "stochastic", steps = steps, solver = obj),
            class = "plant_run")
}

new_plant_run <- function(solver, scenario, record, mutants = NULL) {
  patch <- solver$patch
  structure(list(scenario = scenario,
                 time = solver$time,
                 complete = solver$complete,
                 offspring_production = patch$offspring_production,
                 net_reproduction_ratios = patch$net_reproduction_ratios,
                 event_log = solver$event_log,
                 record = record,
                 mutants = mutants,
                 method = "deterministic",
                 solver = solver),
            class = "plant_run")
}

##' @rdname run
##' @export
run.Individual <- function(x, env = NULL, ..., times = NULL, sizes = NULL,
                           size_name = "height") {
  if (is.null(times) == is.null(sizes)) {
    stop("Give one of `times` and `sizes`", call. = FALSE)
  }
  env <- if (is.null(env)) Environment(x$strategy_name) else live_environment(env)
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
##' @param x A \code{plant_run} run with \code{record = "steps"}, a single
##'   plant's run, or a landscape's run with \code{record = "steps"}.
##' @param ... Unused.
##' @return For a deterministic run, a list: the tidied patch at every step,
##'   \code{steps}, \code{n_spp}, \code{species} and \code{env} (as
##'   \code{\link{tidy_patch}} gives them); then the run's
##'   \code{offspring_production} and \code{net_reproduction_ratios}; \code{p},
##'   the solver's parameters (carrying any refined schedule); \code{events},
##'   the events it was given; \code{event_log}, what each event did; and
##'   \code{control}, the run's \code{Control}, which says how its crowns cast
##'   shade (\code{\link{canopy_profile}} reads it) where the strategies in
##'   \code{p} do not.
##'   For a stochastic run, \code{time}, \code{species} (an array per species,
##'   \code{[variable, time, individual]}, with an \code{is_alive} attribute),
##'   \code{env} and \code{control}. For a single plant, a tibble of its state at each
##'   reported time, as \code{\link{tidy_individual}} gives it. For a
##'   landscape's run, a list of each patch's, by name.
##' @export
tidy.plant_run <- function(x, ...) {
  ## A deterministic result is checked as current first: its steps are the
  ## solver's history, which a later run overwrites.
  if (!identical(x$method, "stochastic")) {
    as_solver(x)
  }
  if (!"steps" %in% x$record) {
    stop("This run kept no steps: run it with record = \"steps\"", call. = FALSE)
  }
  if (identical(x$method, "stochastic")) {
    results <- tidy_stochastic_steps(x$steps)
    results[["control"]] <- x$scenario$control
    return(results)
  }
  solver <- x$solver
  if (length(solver$history) == 0L) {
    stop("This run has taken no step yet: continue it with run(res) first",
         call. = FALSE)
  }
  results <- tidy_patch(lapply(solver$history, "[[", "state"))
  patch <- solver$patch
  results[["offspring_production"]] <- patch$offspring_production
  results[["net_reproduction_ratios"]] <- patch$net_reproduction_ratios
  results[["p"]] <- solver$parameters
  ## Events are given separately from `p`, so `p` alone does not describe the
  ## run: the schedule asked for and what was applied (including how much of
  ## each pulse the soil accepted and shed) travel with the tables.
  results[["events"]] <- solver$events
  results[["event_log"]] <- solver$event_log
  results[["control"]] <- x$scenario$control
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
  solver <- as_solver(resident)
  if (!isTRUE(resident$complete)) {
    stop("The resident has not finished its run; invaders read its whole recording",
         call. = FALSE)
  }
  if (constant) {
    solver$run_mutant_constant(p)
  } else {
    solver$run_mutant(p)
  }
  new_plant_run(solver, resident$scenario, resident$record, mutants = p)
}

## The solver behind a run, for the functions that take either. A result is
## refused once its solver has moved on, since everything read from it would
## describe some later run.
as_solver <- function(x) {
  if (!inherits(x, "plant_run")) {
    return(x)
  }
  if (identical(x$method, "stochastic")) {
    stop("This is a stochastic run, which has no fitness, gradient, export or continuation: ",
         "that needs the deterministic solver", call. = FALSE)
  }
  solver <- x$solver
  time <- tryCatch(solver$time, error = function(e) {
    stop("This result's solver is gone (a result read back with readRDS() keeps ",
         "everything else): ",
         if (is.null(x$mutants)) {
           "run(res, until = res$time) runs its scenario again to where it stood"
         } else {
           "run its resident's scenario again, then run_mutant()"
         }, call. = FALSE)
  })
  if (!identical(time, x$time) ||
      !identical(solver$parameters, x$mutants %||% x$scenario$p)) {
    stop("This result's solver has moved on since it was returned (continued, ",
         "converged or invaded): use the result of whatever moved it", call. = FALSE)
  }
  solver
}
