## Linked patches. A landscape is a set of patches, each a scenario of its own,
## and the flows between them. It is run by the patches' own runners, one
## patch at a time: nothing here knows what a patch holds or what flows.

##' Describe a landscape of linked patches
##'
##' A landscape is a set of patches and the flows that link them. Each patch is
##' a \code{\link{scenario}}: its own species, environment, numerics and events,
##' so patches can differ in anything a single run can. A flow carries something
##' from one patch to another during a run. \code{\link{run}} takes a landscape
##' as it takes a scenario.
##'
##' Nothing here is specific to a model. A landscape knows patches, and flows
##' between them of whatever the patches' environments count as a resource: in
##' TF24 that is water in a soil layer, and \code{resource_flow("ridge",
##' "gully")} passes the part of each rain pulse the ridge's surface cannot
##' hold to the gully.
##'
##' The patches run over one span of time and share its clock, so they must
##' have the same \code{max_patch_lifetime}. They are taken to be of equal
##' area: an amount is per unit ground area in both the patch it leaves and
##' the patch it reaches.
##'
##' Flows must lead one way. A landscape is run patch by patch, each one whole,
##' sources before the patches they feed, so a loop of flows (one patch feeding
##' another that feeds it back) is refused.
##'
##' A landscape is a plain list, and \code{run} checks it again, so a patch or
##' a flow can be changed in place.
##'
##' @param patches A list of \code{\link{scenario}}s, one per patch. Its names
##'   are the patches' names; an unnamed list is named \code{"1"}, \code{"2"},
##'   and so on.
##' @param ... Flows between the patches, as \code{resource_flow} returns
##'   them, or none for patches that run side by side unlinked.
##' @return A list of class \code{plant_landscape} with \code{patches}, the
##'   named scenarios, and \code{flows}, a data frame with one row per flow.
##' @seealso \code{\link{run}}, \code{\link{scenario}},
##'   \code{\link{resource_pulse}}
##' @export
##' @examples
##' p <- default_parameters("TF24", quiet = TRUE)
##' p$max_patch_lifetime <- 2
##' p <- add_strategies(p, trait_matrix(0.0825, "lma"))
##' ## Two patches under the same storms, the upper one on a shallow soil.
##' storms <- events(rainfall_pulse(time = c(0.5, 1, 1.5), depth = 0.06))
##' shallow <- Environment("TF24")
##' shallow$set_soil_layer_widths(c(0.05, 0.1, 0.15, 0.2, 0.25))
##' land <- landscape(
##'   list(ridge = scenario(p, shallow, events = storms),
##'        gully = scenario(p, Environment("TF24"), events = storms)),
##'   resource_flow("ridge", "gully"))
##' \donttest{
##' res <- run(land)
##' res$flows                             # what the ridge passed to the gully
##' res$patches$gully$offspring_production
##' }
landscape <- function(patches, ...) {
  if (is.list(patches) && length(patches) > 0L && is.null(names(patches))) {
    names(patches) <- as.character(seq_along(patches))
  }
  flows <- list(...)
  if (!all(vapply(flows, inherits, logical(1), "plant_flows"))) {
    stop("A landscape's flows are given as resource_flow() returns them",
         call. = FALSE)
  }
  flows <- do.call(rbind, c(list(empty_flows()), flows))
  rownames(flows) <- NULL
  check_landscape(patches, flows)

  structure(list(patches = patches, flows = flows), class = "plant_landscape")
}

## Whether patches and flows make a landscape that can be run. Called where one
## is made and again where it is run, since either part can be changed between.
check_landscape <- function(patches, flows) {
  if (!is.list(patches) || length(patches) == 0L ||
      !all(vapply(patches, inherits, logical(1), "plant_scenario"))) {
    stop("`patches` must be a list of scenarios, as scenario() returns them",
         call. = FALSE)
  }
  if (is.null(names(patches)) || anyNA(names(patches)) ||
      any(names(patches) == "") || anyDuplicated(names(patches)) > 0L) {
    stop("Every patch needs a name of its own", call. = FALSE)
  }
  lifetimes <- vapply(patches, function(sc) sc$p$max_patch_lifetime, numeric(1))
  if (length(unique(lifetimes)) > 1L) {
    stop("The patches of a landscape run over one span of time, and these have ",
         "max_patch_lifetime ", paste(format(lifetimes), collapse = ", "),
         call. = FALSE)
  }

  if (!is.data.frame(flows) ||
      !all(names(empty_flows()) %in% names(flows))) {
    stop("A landscape's flows are given as resource_flow() returns them",
         call. = FALSE)
  }
  check_flows(flows$from, flows$to, flows$resource, flows$into, flows$weight)
  if (!all(flows$flow == "resource")) {
    stop("A flow carries a resource; \"", flows$flow[flows$flow != "resource"][[1L]],
         "\" is not a kind of flow", call. = FALSE)
  }
  unknown <- setdiff(c(flows$from, flows$to), names(patches))
  if (length(unknown) > 0L) {
    stop("A flow names ", paste0("\"", unknown, "\"", collapse = ", "),
         ", which is not a patch of this landscape (",
         paste(names(patches), collapse = ", "), ")", call. = FALSE)
  }
  ## What one pool sheds cannot be passed on more than once over.
  leaving <- tapply(flows$weight, list(flows$from, flows$resource), sum)
  if (any(leaving > 1 + 1e-12, na.rm = TRUE)) {
    over <- which(leaving > 1 + 1e-12, arr.ind = TRUE)[1L, ]
    stop("The flows out of resource ", colnames(leaving)[[over[[2]]]],
         " of patch \"", rownames(leaving)[[over[[1]]]], "\" have weights ",
         "summing to more than 1, which passes on more than the patch sheds",
         call. = FALSE)
  }
  landscape_order(names(patches), flows)
  invisible(TRUE)
}

## What every flow must satisfy, whichever patches it joins.
check_flows <- function(from, to, resource, into, weight) {
  if (!is.character(from) || !is.character(to) || anyNA(from) || anyNA(to)) {
    stop("`from` and `to` are patch names", call. = FALSE)
  }
  if (any(from == to)) {
    stop("A flow links two patches; \"", from[from == to][[1L]],
         "\" cannot flow into itself", call. = FALSE)
  }
  is_index <- function(x) {
    is.numeric(x) && !anyNA(x) && all(x >= 1) && all(x == round(x))
  }
  if (!is_index(resource) || !is_index(into)) {
    stop("`resource` and `into` are indices of resource pools, counted from 1",
         call. = FALSE)
  }
  if (!is.numeric(weight) || anyNA(weight) || any(weight < 0 | weight > 1)) {
    stop("`weight` is a fraction of what is shed, between 0 and 1", call. = FALSE)
  }
}

##' @param from,to The names of the patch the resource leaves and the patch it
##'   reaches.
##' @param resource Index of the resource pool it leaves, in \code{from}'s
##'   environment, as \code{\link{resource_pulse}} counts them: TF24's are its
##'   soil layers, so 1 is the surface.
##' @param into Index of the pool it reaches, in \code{to}'s environment; by
##'   default the same one.
##' @param weight The fraction of what \code{from} sheds that reaches
##'   \code{to}, between 0 and 1. The weights of the flows out of one pool sum
##'   to at most 1; what they leave over leaves the landscape.
##' @details \code{resource_flow} links two patches by what a resource pool
##'   cannot hold. A resource pulse is capped at what its pool can take, and
##'   the rest is shed (\code{\link{resource_pulse}}). With a flow out of that
##'   pool, the shed amount, times \code{weight}, becomes a pulse into the
##'   receiving patch's pool at the same instant, after any pulse the receiver
##'   has of its own at that instant. The receiver may shed it in turn, and
##'   pass it on down its own flows. The arguments are vectorised, so several
##'   flows are one call.
##'
##'   Only what a pulse sheds is passed. A resource a patch loses in any other
##'   way stays lost to the landscape: in TF24, the run-off and drainage of
##'   rain given as a continuous rate (the \code{rainfall} driver) are not
##'   passed on, so a patch whose rain all comes that way passes nothing. Nor
##'   does a flow out of a pool that no pulse reaches, or that the patch's
##'   environment does not have; a flow into a pool the receiving environment
##'   does not have stops the run at that patch.
##' @rdname landscape
##' @export
resource_flow <- function(from, to, resource = 1L, into = resource, weight = 1) {
  if (length(from) == 0L && length(to) == 0L) {
    return(empty_flows())
  }
  n <- max(length(from), length(to), length(resource), length(into),
           length(weight))
  recycle <- function(x, what) {
    if (!(length(x) %in% c(1L, n))) {
      stop("'", what, "' must be length 1 or ", n, call. = FALSE)
    }
    rep_len(x, n)
  }
  from <- recycle(from, "from")
  to <- recycle(to, "to")
  check_flows(from, to, resource, into, weight)
  new_flows(from, to, recycle(as.integer(resource), "resource"),
            recycle(as.integer(into), "into"), recycle(weight, "weight"))
}

## One row per flow. `flow` says what is carried; a resource is the one kind
## so far, and a new kind adds a value here and its own constructor.
new_flows <- function(from, to, resource, into, weight) {
  structure(data.frame(from = from, to = to, flow = rep("resource", length(from)),
                       resource = resource, into = into, weight = weight,
                       stringsAsFactors = FALSE),
            class = c("plant_flows", "data.frame"))
}

empty_flows <- function() {
  new_flows(character(0), character(0), integer(0), integer(0), numeric(0))
}

## The order the patches run in: every patch after the patches that feed it,
## and otherwise as given. Refuses a loop, which no such order satisfies.
landscape_order <- function(patches, flows) {
  order <- character(0)
  remaining <- patches
  while (length(remaining) > 0L) {
    fed <- unique(flows$to[flows$from %in% remaining])
    ready <- setdiff(remaining, fed)
    if (length(ready) == 0L) {
      stop("The flows between ", paste(remaining, collapse = ", "),
           " form a loop. A landscape is run one patch at a time, sources ",
           "first, so its flows must lead one way", call. = FALSE)
    }
    order <- c(order, ready)
    remaining <- setdiff(remaining, ready)
  }
  order
}

##' @export
print.plant_landscape <- function(x, ...) {
  cat(sprintf("<plant_landscape> %d patch(es): %s; %d flow(s)\n",
              length(x$patches), paste(names(x$patches), collapse = ", "),
              nrow(x$flows)))
  invisible(x)
}

##' @rdname run
##' @export
run.plant_landscape <- function(x, ..., until = NULL, from = NULL,
                                record = character(0), refine_schedule = FALSE,
                                method = c("deterministic", "stochastic")) {
  no_more_arguments(...)
  method <- match.arg(method)
  if (!is.null(until) || !is.null(from)) {
    stop("A landscape's patches are each run whole, one after another; it ",
         "cannot take `until` or `from`", call. = FALSE)
  }

  check_landscape(x$patches, x$flows)

  results <- list()
  passed <- list()
  for (patch in landscape_order(names(x$patches), x$flows)) {
    sc <- x$patches[[patch]]
    arriving <- passed[vapply(passed, function(f) f$to[[1L]] == patch, logical(1))]
    if (length(arriving) > 0L) {
      arriving <- do.call(rbind, arriving)
      sc$events <- events(sc$events,
                          resource_pulse(arriving$time, arriving$into,
                                         arriving$amount))
    }
    results[[patch]] <- tryCatch(
      run(sc, record = record, refine_schedule = refine_schedule,
          method = method),
      error = function(e) {
        stop("Patch \"", patch, "\"",
             if (length(arriving) > 0L) ", with what flowed into it",
             ": ", conditionMessage(e), call. = FALSE)
      })
    passed <- c(passed, flows_out_of(patch, results[[patch]], x$flows))
  }

  passed <- if (length(passed) > 0L) do.call(rbind, passed) else {
    data.frame(time = numeric(0), from = character(0), to = character(0),
               resource = integer(0), into = integer(0), amount = numeric(0),
               stringsAsFactors = FALSE)
  }
  passed <- passed[order(passed$time), , drop = FALSE]
  rownames(passed) <- NULL
  structure(list(landscape = x, patches = results[names(x$patches)],
                 flows = passed, method = method),
            class = "plant_landscape_run")
}

## What a patch that has run passes down each of its flows: for every pulse
## into one of its pools that shed something, the flow's share of it. One data
## frame per flow, so each has a single receiver.
flows_out_of <- function(patch, result, flows) {
  flows <- flows[flows$from == patch, , drop = FALSE]
  if (nrow(flows) == 0L) {
    return(list())
  }
  log <- result$event_log
  pulses <- which(log$type == "resource_pulse")
  shed <- vapply(log$applied[pulses], function(a) a[[2L]], numeric(1))
  out <- lapply(seq_len(nrow(flows)), function(i) {
    hit <- log$target_index[pulses] == flows$resource[[i]] & shed > 0
    amount <- flows$weight[[i]] * shed[hit]
    data.frame(time = log$time[pulses][hit], from = rep(patch, sum(hit)),
               to = rep(flows$to[[i]], sum(hit)),
               resource = rep(flows$resource[[i]], sum(hit)),
               into = rep(flows$into[[i]], sum(hit)), amount = amount,
               stringsAsFactors = FALSE)[amount > 0, , drop = FALSE]
  })
  out[vapply(out, nrow, integer(1)) > 0L]
}

##' @export
print.plant_landscape_run <- function(x, ...) {
  cat(sprintf("<plant_landscape_run> %d patch(es)%s, %d transfer(s) between them\n",
              length(x$patches),
              if (identical(x$method, "stochastic")) " (stochastic)" else "",
              nrow(x$flows)))
  for (patch in names(x$patches)) {
    cat("  ", patch, ": ", sep = "")
    print(x$patches[[patch]])
  }
  invisible(x)
}

##' @rdname tidy.plant_run
##' @export
tidy.plant_landscape_run <- function(x, ...) {
  lapply(x$patches, tidy, ...)
}
