##' Discrete events applied during an SCM run.
##'
##' An event is a \code{(time, action)} pair applied between solver legs: the
##' integrator is stopped at the event time, the action changes the patch, the
##' system is recomputed, and integration resumes from the new state. Node
##' introduction has always worked this way; these functions make the same
##' mechanism available for rainfall pulses, thinning and heat damage.
##'
##' Each event carries when it happens, what kind of thing it is, what it acts
##' on, and the values it needs. The target is \code{"environment"} (the
##' abiotic state), \code{"patch"} (every species) or \code{"species"} (one,
##' named by \code{species}). There is deliberately no per-cohort target: a
##' cohort has no stable address across a run, because nodes are appended and
##' never removed and schedule refinement changes how many exist. Selecting
##' particular cohorts is expressed as a size range instead, which is both well
##' defined and what thinning actually needs.
##'
##' Events are instantaneous \emph{to the solver}: patch time does not advance
##' across one. An action is still free to reach its answer by integrating its
##' own fast sub-model over a nominal duration — \code{\link{heat_damage}} does
##' — with the patch's demography frozen. The solver sees one jump either way.
##'
##' Events sharing a time are applied in a fixed order: environment events
##' first, then demographic ones, then node introductions — so a node
##' introduced at that instant starts life in the post-event environment.
##'
##' Because an event is also a stop time for the integrator, adding one changes
##' the adaptive step sequence. A run with events legitimately differs from one
##' without, at solver tolerance, even away from the events themselves.
##'
##' What each event actually did — as against what was asked of it — is
##' recorded, and readable afterwards as \code{scm$event_log}. The two differ
##' routinely: a pulse is capped at what the soil can hold, and thinning a size
##' class removes whatever was in it.
##'
##' @param ... For \code{events}, objects returned by the individual event
##'   constructors, or whole \code{Events} objects; these are concatenated. Each
##'   constructor is vectorised over its arguments, so a whole series of
##'   rainfall pulses is one call.
##' @return An \code{Events} object: a list with \code{time}, \code{type},
##'   \code{target}, \code{target_index} and \code{params}, in schedule order.
##' @rdname events
##' @export
##' @examples
##' p <- scm_base_parameters("FF16")
##' p <- add_strategies(p, trait_matrix(1, "lma"))
##' ev <- events(
##'   events_default(p),
##'   thinning(time = 20, fraction = 0.3)
##' )
events <- function(...) {
  parts <- list(...)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) {
    return(empty_events())
  }
  ev <- do.call(join_event_rows, parts)
  ## Sorting here is a convenience for anyone reading the object; the C++
  ## queue re-sorts on insertion and is the authority on ties.
  i <- order(ev$time)
  Events(time = as.numeric(ev$time[i]),
         type = as.character(ev$type[i]),
         target = as.character(ev$target[i]),
         target_index = as.integer(ev$target_index[i]),
         params = ev$params[i])
}

## The "no events supplied" signal: the C++ side falls back to
## p$node_schedule_times when the list is empty.
empty_events <- function() {
  Events(time = numeric(0), type = character(0), target = character(0),
         target_index = integer(0), params = list())
}

## Concatenate event rows. unlist() on an empty list gives NULL, which would
## silently turn an empty schedule into a malformed one, so the empty case is
## spelled out with the right types.
join_event_rows <- function(...) {
  parts <- list(...)
  pull <- function(field, empty) {
    if (length(parts) == 0L) return(empty)
    out <- unlist(lapply(parts, `[[`, field), use.names = FALSE)
    if (is.null(out)) empty else out
  }
  list(time = pull("time", numeric(0)),
       type = pull("type", character(0)),
       target = pull("target", character(0)),
       target_index = pull("target_index", integer(0)),
       params = if (length(parts) == 0L) list() else
         unlist(lapply(parts, `[[`, "params"), recursive = FALSE))
}

## Build one event type's rows. Vectorised over time and over each parameter;
## all must be length 1 or a common length n.
event_rows <- function(type, time, target, target_index = 1L, params = list()) {
  n <- max(c(length(time), vapply(params, length, integer(1)), 1L))
  recycle <- function(x, what) {
    if (length(x) == n) return(x)
    if (length(x) == 1L) return(rep(x, n))
    stop(sprintf("'%s' must be length 1 or %d, not %d", what, n, length(x)),
         call. = FALSE)
  }
  time <- recycle(time, "time")
  target <- recycle(target, "target")
  target_index <- recycle(target_index, "species")
  params <- lapply(seq_along(params),
                   function(j) recycle(params[[j]], names(params)[[j]]))
  list(time = as.numeric(time),
       type = rep(type, n),
       target = as.character(target),
       target_index = as.integer(target_index),
       ## One numeric vector per event, in the order the C++ action reads them.
       params = lapply(seq_len(n), function(i) {
         as.numeric(vapply(params, function(p) p[[i]], numeric(1)))
       }))
}

## "patch" unless a species was named, in which case "species".
scope_of <- function(species) {
  if (is.null(species)) "patch" else "species"
}

##' @details \code{events_default(p)} is the schedule a run gets when no events
##'   are supplied: the node introductions from \code{p$node_schedule_times} and
##'   nothing else. Start from it when adding events to an otherwise ordinary
##'   run — \code{events(events_default(p), rainfall_pulse(...))} — or pass it
##'   on its own, which reproduces the default run exactly.
##' @rdname events
##' @export
events_default <- function(p) {
  events(node_introductions(p))
}

##' @param p Parameters object, whose \code{node_schedule_times} supply the
##'   introduction times.
##' @rdname events
##' @export
node_introductions <- function(p) {
  times <- p$node_schedule_times
  parts <- lapply(seq_along(times), function(i) {
    event_rows("node_introduction", time = times[[i]],
               target = "species", target_index = i)
  })
  do.call(join_event_rows, parts)
}

##' @param time Event time(s), in years of patch age.
##' @param depth Depth of water delivered, in m (so 13 mm is \code{0.013}).
##'   What the surface layer cannot hold is recorded as runoff rather than
##'   forced in; see \code{scm$event_log}.
##' @rdname events
##' @export
rainfall_pulse <- function(time, depth) {
  event_rows("rainfall_pulse", time = time, target = "environment",
             params = list(depth = depth))
}

##' @param fraction Fraction of individuals removed, in \code{[0, 1)}.
##' @param height_min,height_max Only individuals whose height falls in this
##'   band are removed (m). The defaults take the whole stand.
##' @param species Index of the species to act on; \code{NULL} (the default)
##'   acts on every species in the patch.
##' @rdname events
##' @export
thinning <- function(time, fraction, height_min = 0, height_max = Inf,
                     species = NULL) {
  event_rows("thinning", time = time, target = scope_of(species),
             target_index = if (is.null(species)) 1L else species,
             params = list(fraction = fraction, height_min = height_min,
                           height_max = height_max))
}

##' @details \code{harvest} and \code{partial_disturbance} are \code{thinning}
##'   under two names that read better at their own call sites: a harvest takes
##'   everything above a size, a partial disturbance takes a fraction of
##'   everything. They are the same action — the only difference is which
##'   individuals are selected.
##' @rdname events
##' @export
harvest <- function(time, fraction, height_min = 0, species = NULL) {
  thinning(time, fraction = fraction, height_min = height_min,
           species = species)
}

##' @rdname events
##' @export
partial_disturbance <- function(time, fraction, species = NULL) {
  thinning(time, fraction = fraction, species = species)
}

##' @param temperature Peak air temperature reached during the event (deg C).
##' @param duration Nominal duration of the event (years; a fortnight is
##'   \code{14 / 365}). The action sub-integrates over it at half-hourly steps;
##'   patch time does not advance.
##' @param temperature_crit Temperature above which damage accrues (deg C).
##' @param sensitivity Damage accrued per degree-year above
##'   \code{temperature_crit}.
##' @rdname events
##' @export
heat_damage <- function(time, temperature, duration = 14 / 365,
                        temperature_crit = 40, sensitivity = 1,
                        species = NULL) {
  event_rows("heat_damage", time = time, target = scope_of(species),
             target_index = if (is.null(species)) 1L else species,
             params = list(temperature = temperature, duration = duration,
                           temperature_crit = temperature_crit,
                           sensitivity = sensitivity))
}
