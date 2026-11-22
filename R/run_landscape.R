## Linked patches. A landscape is a set of patches, each a scenario of its own,
## and the flows between them. It is run by the patches' own runners, one
## patch at a time. The landscape knows patches and a table of flows; what a
## kind of flow carries, and how, is that kind's own (see flow_kinds()).

##' Describe a landscape of linked patches
##'
##' A landscape is a set of patches and the flows that link them. Each patch is
##' a \code{\link{scenario}}: its own species, environment, numerics and events,
##' so patches can differ in anything a single run can. A flow carries something
##' from one patch to another during a run. \code{\link{run}} takes a landscape
##' as it takes a scenario.
##'
##' Nothing here is specific to a model, or to what flows. Every flow is one
##' row of the same table: the patch it leaves (\code{from}) and the patch it
##' reaches (\code{to}); its \code{type}, which says what it carries; which
##' part of the source it comes \code{out_of} and which part of the receiver
##' it goes \code{into}, each an index (of a resource pool of the environment,
##' or of a species, as the type reads them); and its \code{weight}, the share
##' of what leaves that reaches \code{to}. A kind of flow is a function that
##' writes rows of this table, so a landscape can mix kinds.
##'
##' One kind exists so far. \code{resource_flow} carries whatever the patches'
##' environments count as a resource: in TF24 that is water in a soil layer,
##' and \code{resource_flow("ridge", "gully")} passes the part of each rain
##' pulse the ridge's surface cannot hold to the gully.
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
##'   named scenarios, and \code{flows}, a data frame with one row per flow:
##'   \code{from}, \code{to}, \code{type}, \code{out_of}, \code{into} and
##'   \code{weight}.
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
    stop("A landscape's flows are rows of its flows table, as resource_flow() ",
         "returns them", call. = FALSE)
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
    stop("A landscape's flows are rows of its flows table, as resource_flow() ",
         "returns them", call. = FALSE)
  }
  check_flows(flows)
  unknown <- setdiff(c(flows$from, flows$to), names(patches))
  if (length(unknown) > 0L) {
    stop("A flow names ", paste0("\"", unknown, "\"", collapse = ", "),
         ", which is not a patch of this landscape (",
         paste(names(patches), collapse = ", "), ")", call. = FALSE)
  }
  landscape_order(names(patches), flows)
  invisible(TRUE)
}

## What the rows of a flows table must satisfy, whichever patches they join:
## what holds of any flow, then what each kind asks of its own rows.
check_flows <- function(flows) {
  if (nrow(flows) == 0L) {
    return(invisible(TRUE))
  }
  if (!is.character(flows$from) || !is.character(flows$to) ||
      anyNA(flows$from) || anyNA(flows$to)) {
    stop("A flow's `from` and `to` are patch names", call. = FALSE)
  }
  kinds <- flow_kinds()
  unknown <- setdiff(flows$type, names(kinds))
  if (length(unknown) > 0L) {
    stop("\"", unknown[[1L]], "\" is not a kind of flow; the kinds are ",
         paste(names(kinds), collapse = ", "), call. = FALSE)
  }
  is_index <- function(x) {
    is.numeric(x) && !anyNA(x) && all(x >= 1) && all(x == round(x))
  }
  if (!is_index(flows$out_of) || !is_index(flows$into)) {
    stop("A flow's `out_of` and `into` are indices, counted from 1",
         call. = FALSE)
  }
  if (!is.numeric(flows$weight) || anyNA(flows$weight) ||
      any(flows$weight < 0 | flows$weight > 1)) {
    stop("A flow's `weight` is a share of what leaves, between 0 and 1",
         call. = FALSE)
  }
  ## What leaves one part of a patch cannot be passed on more than once over.
  source <- paste(flows$type, flows$out_of, "of patch", flows$from)
  leaving <- tapply(flows$weight, source, sum)
  if (any(leaving > 1 + 1e-12)) {
    stop("The flows out of ", names(leaving)[leaving > 1 + 1e-12][[1L]],
         " have weights summing to more than 1, which passes on more than ",
         "leaves it", call. = FALSE)
  }
  for (kind in unique(flows$type)) {
    kinds[[kind]]$check(flows[flows$type == kind, , drop = FALSE])
  }
  invisible(TRUE)
}

## The kinds of flow, by the name a flow's `type` holds. A kind is three
## functions, and adding one is adding an entry here and a constructor of its
## rows:
##   check(flows)             -- refuses rows of this kind that make no sense
##                               for it, beyond what check_flows() asks of all;
##   leaving(result, out_of)  -- what a patch that has run let go from the part
##                               `out_of` indexes: a data frame of `time` and
##                               `amount`, before any flow's weight;
##   arriving(scenario, arrivals) -- the receiving patch's scenario with the
##                               amounts in `arrivals` (`time`, `into`,
##                               `amount`) added to it, before it runs.
## The landscape's own run reads nothing else of a kind, so what one carries
## and how it enters a patch stay here.
flow_kinds <- function() {
  list(resource = list(check = resource_flow_check,
                       leaving = resource_leaving,
                       arriving = resource_arriving))
}

## A resource leaves one patch's environment for another's.
resource_flow_check <- function(flows) {
  if (any(flows$from == flows$to)) {
    stop("A resource flows between two patches; \"",
         flows$from[flows$from == flows$to][[1L]],
         "\" cannot flow into itself", call. = FALSE)
  }
}

## What each pulse into the pool shed: the amount its pool could not hold.
resource_leaving <- function(result, out_of) {
  log <- result$event_log
  pulses <- which(log$type == "resource_pulse" & log$target_index == out_of)
  shed <- vapply(log$applied[pulses], function(a) a[[2L]], numeric(1))
  data.frame(time = log$time[pulses], amount = shed)
}

## Each amount arrives as a pulse into the receiver's pool at its instant,
## after the receiver's own events there.
resource_arriving <- function(scenario, arrivals) {
  scenario$events <- events(scenario$events,
                            resource_pulse(arrivals$time, arrivals$into,
                                           arrivals$amount))
  scenario
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
  flows <- new_flows(recycle(from, "from"), recycle(to, "to"), "resource",
                     recycle(resource, "resource"), recycle(into, "into"),
                     recycle(weight, "weight"))
  check_flows(flows)
  flows$out_of <- as.integer(flows$out_of)
  flows$into <- as.integer(flows$into)
  flows
}

## The flows table: one row per flow, the same columns whatever is carried.
new_flows <- function(from, to, type, out_of, into, weight) {
  structure(data.frame(from = from, to = to, type = rep(type, length(from)),
                       out_of = out_of, into = into, weight = weight,
                       stringsAsFactors = FALSE),
            class = c("plant_flows", "data.frame"))
}

empty_flows <- function() {
  new_flows(character(0), character(0), character(0), integer(0), integer(0),
            numeric(0))
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

  kinds <- flow_kinds()
  results <- list()
  passed <- empty_transfers()
  for (patch in landscape_order(names(x$patches), x$flows)) {
    sc <- x$patches[[patch]]
    arriving <- passed[passed$to == patch, , drop = FALSE]
    for (kind in unique(arriving$type)) {
      sc <- kinds[[kind]]$arriving(sc, arriving[arriving$type == kind, , drop = FALSE])
    }
    results[[patch]] <- tryCatch(
      run(sc, record = record, refine_schedule = refine_schedule,
          method = method),
      error = function(e) {
        stop("Patch \"", patch, "\"",
             if (nrow(arriving) > 0L) ", with what flowed into it",
             ": ", conditionMessage(e), call. = FALSE)
      })
    passed <- rbind(passed, transfers_out_of(patch, results[[patch]], x$flows, kinds))
  }

  passed <- passed[order(passed$time), , drop = FALSE]
  rownames(passed) <- NULL
  structure(list(landscape = x, patches = results[names(x$patches)],
                 flows = passed, method = method),
            class = "plant_landscape_run")
}

## Every amount passed between patches: when, from where to where, and how
## much reached the receiver.
empty_transfers <- function() {
  data.frame(time = numeric(0), from = character(0), to = character(0),
             type = character(0), out_of = integer(0), into = integer(0),
             amount = numeric(0), stringsAsFactors = FALSE)
}

## What a patch that has run passes down each of its flows: the flow's share
## of each amount its kind says left the part it comes out of.
transfers_out_of <- function(patch, result, flows, kinds) {
  flows <- flows[flows$from == patch, , drop = FALSE]
  out <- lapply(seq_len(nrow(flows)), function(i) {
    left <- kinds[[flows$type[[i]]]]$leaving(result, flows$out_of[[i]])
    if (!all(is.finite(left$amount)) || any(left$amount < 0)) {
      stop("Patch \"", patch, "\" reported an amount leaving it by a ",
           flows$type[[i]], " flow that is not a finite, non-negative number",
           call. = FALSE)
    }
    amount <- flows$weight[[i]] * left$amount
    n <- nrow(left)
    data.frame(time = left$time, from = rep(patch, n),
               to = rep(flows$to[[i]], n), type = rep(flows$type[[i]], n),
               out_of = rep(flows$out_of[[i]], n),
               into = rep(flows$into[[i]], n), amount = amount,
               stringsAsFactors = FALSE)[amount > 0, , drop = FALSE]
  })
  do.call(rbind, c(list(empty_transfers()), out))
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
