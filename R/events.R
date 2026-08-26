##' Discrete events applied during an SCM run.
##'
##' An event is a \code{(time, action)} pair applied between solver legs: the
##' integrator is stopped at the event time, the action changes the patch, and
##' integration resumes from the new state. Node introduction has always worked
##' this way; these functions make the same mechanism available for rainfall
##' pulses, harvest, partial disturbance and temperature extremes.
##'
##' Events are instantaneous \emph{to the solver} — patch time does not advance
##' across one. An action is still free to reach its answer by integrating its
##' own fast sub-model over a nominal duration (a fortnight of heat, say) with
##' the patch's demography frozen; that is a private matter for the action, and
##' the solver sees one jump either way.
##'
##' Events sharing a time are applied in a fixed order: environment events
##' (rainfall, temperature) first, then demographic removals (harvest, partial
##' disturbance), then node introductions — so a node introduced at that
##' instant starts life in the post-event environment.
##'
##' @param ... For \code{events}, objects returned by the individual event
##'   constructors; these are concatenated. Each constructor is vectorised over
##'   its arguments, so a whole series of rainfall pulses is one call.
##' @return An \code{Events} object: a list with \code{time}, \code{type},
##'   \code{species_index} and \code{params}, in schedule order.
##' @rdname events
##' @export
##' @examples
##' p <- scm_base_parameters("FF16")
##' p <- add_strategies(p, trait_matrix(1, "lma"))
##' ev <- events(
##'   node_introductions(p),
##'   partial_disturbance(time = 20, fraction = 0.3)
##' )
events <- function(...) {
  parts <- list(...)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) {
    return(empty_events())
  }
  ev <- list(
    time          = unlist(lapply(parts, function(x) x$time), use.names = FALSE),
    type          = unlist(lapply(parts, function(x) x$type), use.names = FALSE),
    species_index = unlist(lapply(parts, function(x) x$species_index),
                           use.names = FALSE),
    params        = unlist(lapply(parts, function(x) x$params), recursive = FALSE)
  )
  ## Sorting here is a convenience for anyone reading the object; the C++
  ## queue re-sorts on insertion and is the authority on ties.
  i <- order(ev$time)
  Events(time = as.numeric(ev$time[i]),
         type = as.character(ev$type[i]),
         species_index = as.integer(ev$species_index[i]),
         params = ev$params[i])
}

## The "no events supplied" signal: the C++ side falls back to
## p$node_schedule_times when the list is empty.
empty_events <- function() {
  Events(time = numeric(0), type = character(0),
         species_index = integer(0), params = list())
}

## Build one event type's rows. Vectorised over time and over each parameter;
## all must be length 1 or a common length n.
event_rows <- function(type, time, species_index = 1L, params = list()) {
  n <- max(c(length(time), vapply(params, length, integer(1)), 1L))
  recycle <- function(x, what) {
    if (length(x) == n) return(x)
    if (length(x) == 1L) return(rep(x, n))
    stop(sprintf("'%s' must be length 1 or %d, not %d", what, n, length(x)),
         call. = FALSE)
  }
  time <- recycle(time, "time")
  species_index <- recycle(species_index, "species_index")
  params <- lapply(seq_along(params), function(j) recycle(params[[j]], names(params)[[j]]))
  list(time = as.numeric(time),
       type = rep(type, n),
       species_index = as.integer(species_index),
       ## One numeric vector per event, in the order the C++ action reads them.
       params = lapply(seq_len(n), function(i) {
         as.numeric(vapply(params, function(p) p[[i]], numeric(1)))
       }))
}

##' @param p Parameters object, whose \code{node_schedule_times} supply the
##'   introduction times.
##' @rdname events
##' @export
node_introductions <- function(p) {
  times <- p$node_schedule_times
  parts <- lapply(seq_along(times), function(i) {
    event_rows("node_introduction", time = times[[i]], species_index = i)
  })
  do.call(join_event_rows, parts)
}

join_event_rows <- function(...) {
  parts <- list(...)
  list(time = unlist(lapply(parts, `[[`, "time"), use.names = FALSE),
       type = unlist(lapply(parts, `[[`, "type"), use.names = FALSE),
       species_index = unlist(lapply(parts, `[[`, "species_index"), use.names = FALSE),
       params = unlist(lapply(parts, `[[`, "params"), recursive = FALSE))
}

##' @param time Event time(s), in years of patch age.
##' @param depth Depth of water delivered, in m (so 13 mm is \code{0.013}).
##'   The amount layer 0 can accept is capped at its free capacity; the excess
##'   is recorded as pulse runoff.
##' @rdname events
##' @export
rainfall_pulse <- function(time, depth) {
  event_rows("rainfall_pulse", time = time, params = list(depth = depth))
}

##' @param fraction Fraction of individuals removed, in \code{[0, 1)}.
##' @param height_min Only individuals taller than this are harvested (m).
##' @rdname events
##' @export
harvest <- function(time, fraction, height_min = 0) {
  event_rows("harvest", time = time,
             params = list(fraction = fraction, height_min = height_min))
}

##' @rdname events
##' @export
partial_disturbance <- function(time, fraction) {
  event_rows("partial_disturbance", time = time,
             params = list(fraction = fraction))
}

##' @param temperature Peak air temperature reached during the event (deg C).
##' @param duration Nominal duration of the event (years; a fortnight is
##'   \code{14 / 365}). Used only by the action's internal sub-model — patch
##'   time does not advance across the event.
##' @param temperature_crit Temperature above which damage accrues (deg C).
##' @param sensitivity Damage accrued per degree-year above
##'   \code{temperature_crit}.
##' @rdname events
##' @export
temperature_extreme <- function(time, temperature, duration = 14 / 365,
                                temperature_crit = 40, sensitivity = 1) {
  event_rows("temperature_extreme", time = time,
             params = list(temperature = temperature, duration = duration,
                           temperature_crit = temperature_crit,
                           sensitivity = sensitivity))
}
