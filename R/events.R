##' Discrete events applied during an SCM run.
##'
##' An event is a \code{(time, action)} pair applied between solver legs: the
##' integrator is stopped at the event time, the action changes the patch, the
##' system is recomputed, and integration resumes from the new state. Node
##' introduction has always worked this way; these functions make the same
##' mechanism available for resource pulses, harvest and climate extremes.
##'
##' The vocabulary here is deliberately model-agnostic, because the machinery is
##' shared by every strategy and environment and nothing about it is specific to
##' plants: a resource pulse is water in TF24 and could be anything countable in
##' a size-structured animal model, and a climate extreme is heat in one model
##' and could be cold or salinity in another. Names that are only true of one
##' model live with that model — \code{\link{rainfall_pulse}} is a resource
##' pulse of water into TF24's surface soil layer.
##'
##' Each event carries when it happens, what kind of thing it is, what it acts
##' on, and the values it needs. The target is \code{"environment"} (its
##' resource pools and drivers), \code{"patch"} (every species) or
##' \code{"species"} (one, named by \code{species}). There is deliberately no
##' per-cohort target: a cohort has no stable address across a run, because
##' nodes are appended and never removed and schedule refinement changes how
##' many exist. Selecting particular cohorts is expressed as a size range
##' instead, which is both well defined and what size-selective removal needs.
##'
##' Events are instantaneous \emph{to the solver}: patch time does not advance
##' across one. An action is still free to reach its answer by integrating its
##' own fast sub-model over a nominal duration — \code{\link{climate_extreme}}
##' does — with demography frozen. The solver sees one jump either way.
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
##' routinely: a pulse is capped at what the pool can hold, and harvesting a
##' size class removes whatever was in it.
##'
##' @param ... For \code{events}, objects returned by the individual event
##'   constructors, or whole \code{Events} objects; these are concatenated. Each
##'   constructor is vectorised over its arguments, so a whole series of
##'   pulses is one call.
##' @return An \code{Events} object: a list with \code{time}, \code{type},
##'   \code{target}, \code{target_index} and \code{params}, in schedule order.
##' @rdname events
##' @export
##' @examples
##' p <- scm_base_parameters("FF16")
##' p <- add_strategies(p, trait_matrix(1, "lma"))
##' ev <- events(
##'   events_default(p),
##'   harvest(time = 20, fraction = 0.3)
##' )
events <- function(...) {
  parts <- list(...)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) {
    return(empty_events())
  }
  ev <- do.call(join_event_rows, parts)
  ## Sort the way the queue applies them -- by time, then by type -- so the
  ## object you read is in the order things actually happen. Sorting by time
  ## alone would list a different order from the one that runs, which is the
  ## kind of discrepancy nobody checks until it matters. `order` is stable, so
  ## events of one type at one instant keep the order they were given in, and
  ## that matters: two pulses at one time are capped in sequence against the
  ## same pool, so the order decides which is credited with the water.
  i <- order(ev$time, match(ev$type, event_type_order()))
  Events(time = as.numeric(ev$time[i]),
         type = as.character(ev$type[i]),
         target = as.character(ev$target[i]),
         target_index = as.integer(ev$target_index[i]),
         params = ev$params[i])
}

## The within-time application order, mirroring EventType in node_schedule.h:
## environment first (it sets the conditions), then removals, then node
## introduction last so a newborn starts in the post-event environment.
## Kept in step with the C++ enum by the round-trip test in test-events.R.
event_type_order <- function() {
  c("resource_pulse", "climate_extreme", "harvest", "node_introduction")
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
##'   run — \code{events(events_default(p), harvest(...))} — or pass it
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
##' @param resource Index of the environment's resource pool to add to. What a
##'   resource is depends on the environment: TF24's are its soil layers, so
##'   \code{resource = 1} is the surface. See \code{\link{rainfall_pulse}} for
##'   that case under a name that reads correctly.
##' @param amount How much to add, in whatever the environment measures that
##'   resource in. What the pool cannot hold is shed rather than forced in;
##'   \code{scm$event_log} reports how much of each.
##' @rdname events
##' @export
resource_pulse <- function(time, resource, amount) {
  event_rows("resource_pulse", time = time, target = "environment",
             target_index = resource, params = list(amount = amount))
}

##' @param fraction Fraction of individuals removed, in \code{[0, 1)}.
##' @param size_min,size_max Only individuals whose size falls in this band are
##'   removed. Size is the individual's size coordinate — height, for plant's
##'   models. The defaults take the whole population.
##' @param species Index of the species to act on; \code{NULL} (the default)
##'   acts on every species in the patch.
##' @details \code{harvest} covers removal generally, not only the forestry
##'   sense: leave the size band at its defaults and it is an across-the-board
##'   knock-down; set \code{size_min} and it takes everything above a size; set
##'   both and it thins one size class. They are one action, differing only in
##'   which individuals are selected.
##' @rdname events
##' @export
harvest <- function(time, fraction, size_min = 0, size_max = Inf,
                    species = NULL) {
  event_rows("harvest", time = time, target = scope_of(species),
             target_index = if (is.null(species)) 1L else species,
             params = list(fraction = fraction, size_min = size_min,
                           size_max = size_max))
}

##' @param intensity Peak level reached during the episode, in whatever the
##'   model measures the stressor in (deg C for a heatwave, say).
##' @param duration Nominal duration of the episode (years; a fortnight is
##'   \code{14 / 365}). The action sub-integrates over it at half-hourly steps;
##'   patch time does not advance.
##' @param threshold Level above which dose accrues; below it the episode does
##'   nothing.
##' @param sensitivity Dose accrued per unit-year above \code{threshold}.
##' @rdname events
##' @export
climate_extreme <- function(time, intensity, duration = 14 / 365,
                            threshold = 40, sensitivity = 1,
                            species = NULL) {
  event_rows("climate_extreme", time = time, target = scope_of(species),
             target_index = if (is.null(species)) 1L else species,
             params = list(intensity = intensity, duration = duration,
                           threshold = threshold,
                           sensitivity = sensitivity))
}
