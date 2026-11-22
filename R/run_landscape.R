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
##' it goes \code{into}, each an index that the type reads; and its
##' \code{weight}, how strongly it joins the two.
##'
##' \code{flow} writes rows of that table. One kind of flow exists so far:
##' \code{"resource"}, which carries whatever the patches' environments count
##' as a resource. In TF24 that is water in a soil layer, and
##' \code{flow("ridge", "gully", "resource")} passes the part of each rain
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
##' @param ... For \code{landscape}, the flows between the patches, as
##'   \code{flow} returns them, or none for patches that run side by side
##'   unlinked. For \code{flow}, further values that its kind takes; a
##'   resource flow takes none.
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
##'   flow("ridge", "gully", "resource"))
##' \donttest{
##' res <- run(land)
##' res$transfers                         # what the ridge passed to the gully
##' res$patches$gully$offspring_production
##' }
landscape <- function(patches, ...) {
  if (is.list(patches) && length(patches) > 0L && is.null(names(patches))) {
    names(patches) <- as.character(seq_along(patches))
  }
  flows <- list(...)
  if (!all(vapply(flows, inherits, logical(1), "plant_flows"))) {
    stop("A landscape's flows are rows of its flows table, as flow() returns ",
         "them", call. = FALSE)
  }
  flows <- bind_flows(c(list(empty_flows()), flows))
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
    stop("A landscape's flows are rows of its flows table, as flow() returns ",
         "them", call. = FALSE)
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
  if (!is.character(flows$type) || anyNA(flows$type)) {
    stop("A flow's `type` is the name of its kind: one of ",
         paste(names(kinds), collapse = ", "), call. = FALSE)
  }
  unknown <- setdiff(flows$type, names(kinds))
  if (length(unknown) > 0L) {
    stop("\"", unknown[[1L]], "\" is not a kind of flow; the kinds are ",
         paste(names(kinds), collapse = ", "), call. = FALSE)
  }
  if (!is_index(flows$out_of) || !is_index(flows$into)) {
    stop("A flow's `out_of` and `into` are indices, counted from 1",
         call. = FALSE)
  }
  if (!is.numeric(flows$weight) || !all(is.finite(flows$weight)) ||
      any(flows$weight < 0)) {
    stop("A flow's `weight` is a finite number, zero or more", call. = FALSE)
  }
  for (kind in unique(flows$type)) {
    kinds[[kind]]$check(flows[flows$type == kind, , drop = FALSE])
  }
  invisible(TRUE)
}

## A whole number an index can be: finite, at least 1 and within an integer.
is_index <- function(x) {
  is.numeric(x) && all(is.finite(x)) && all(x >= 1) && all(x == round(x)) &&
    all(x <= .Machine$integer.max)
}

## Flows tables stacked into one. A kind may give its rows columns of its
## own, which the other kinds' rows hold as NA.
bind_flows <- function(tables) {
  columns <- unique(unlist(lapply(tables, names)))
  tables <- lapply(tables, function(x) {
    x <- as.data.frame(x)
    for (column in setdiff(columns, names(x))) {
      x[[column]] <- rep(NA, nrow(x))
    }
    x[columns]
  })
  flows <- do.call(rbind, tables)
  rownames(flows) <- NULL
  structure(flows, class = c("plant_flows", "data.frame"))
}

## The kinds of flow, by the name a flow's `type` holds. A kind is three
## functions, and the landscape's run reads a kind through these alone, with
## `columns`, the names of any columns its rows carry beyond the table's own
## (flow() takes a value for each and refuses any other):
##   check(flows)            -- refuses rows of this kind that make no sense
##                              for it, beyond what check_flows() asks of all;
##   leaving(result, flow)   -- what a patch that has run let go by `flow`, one
##                              row of the flows table: a data frame of `time`
##                              and `amount`, before the flow's weight;
##   arriving(scenario, arrivals) -- the receiving patch's scenario with the
##                              amounts in `arrivals` (rows of the transfers
##                              table, in time order) added to it, before it
##                              runs.
## A new kind is an entry here; flow() writes its rows. What the run
## cannot do for any kind is in run.plant_landscape(): each patch is run whole,
## sources first, so a flow that feeds back, into its own patch included, is
## refused there whatever its kind.
flow_kinds <- function() {
  list(resource = list(columns = character(0),
                       check = resource_flow_check,
                       leaving = resource_leaving,
                       arriving = resource_arriving))
}

## A resource leaves one patch's environment for another's, and no more of it
## arrives than left: a weight is a share, and the shares out of one pool sum
## to at most 1.
resource_flow_check <- function(flows) {
  if (any(flows$from == flows$to)) {
    stop("A resource flows between two patches; \"",
         flows$from[flows$from == flows$to][[1L]],
         "\" cannot flow into itself", call. = FALSE)
  }
  if (any(flows$weight > 1)) {
    stop("A resource flow's `weight` is a share of what is shed, between 0 ",
         "and 1", call. = FALSE)
  }
  source <- paste("resource", flows$out_of, "of patch", paste0("\"", flows$from, "\""))
  leaving <- tapply(flows$weight, source, sum)
  if (any(leaving > 1 + 1e-12)) {
    stop("The flows out of ", names(leaving)[leaving > 1 + 1e-12][[1L]],
         " have weights summing to more than 1, which passes on more than ",
         "the patch sheds", call. = FALSE)
  }
}

## What each pulse into the pool shed: the amount its pool could not hold, at
## the pulses that shed any.
resource_leaving <- function(result, flow) {
  log <- result$event_log
  pulses <- which(log$type == "resource_pulse" &
                    log$target_index == flow$out_of)
  shed <- vapply(log$applied[pulses], function(a) a[[2L]], numeric(1))
  data.frame(time = log$time[pulses], amount = shed)[shed > 0, , drop = FALSE]
}

## Each amount arrives as a pulse into the receiver's pool at its instant,
## after the receiver's own events there.
resource_arriving <- function(scenario, arrivals) {
  scenario$events <- events(scenario$events,
                            resource_pulse(arrivals$time, arrivals$into,
                                           arrivals$amount))
  scenario
}

##' @param from,to The names of the patch the flow leaves and the patch it
##'   reaches.
##' @param type The kind of flow, which says what it carries:
##'   \code{"resource"}.
##' @param out_of Which part of \code{from}'s patch the flow comes out of, an
##'   index that its kind reads. For a resource flow, a resource pool of the
##'   environment, as \code{\link{resource_pulse}} counts them: TF24's are its
##'   soil layers, so 1 is the surface.
##' @param into Which part of \code{to}'s patch it goes into, read the same
##'   way; by default the same index.
##' @param weight How strongly the flow joins the two patches: what leaves
##'   \code{from} is multiplied by it on the way to \code{to}. For a resource
##'   flow, the share of what \code{from} sheds that reaches \code{to},
##'   between 0 and 1.
##' @details \code{flow}'s arguments are vectorised over everything but
##'   \code{type}, so several flows of one kind are one call.
##' @section Resource flows:
##'   A resource flow links two patches by what a resource pool cannot hold. A
##'   resource pulse is capped at what its pool can take, and the rest is shed
##'   (\code{\link{resource_pulse}}). With a flow out of that pool, the shed
##'   amount, times \code{weight}, becomes a pulse into the receiving patch's
##'   pool at the same instant, after any pulse the receiver has of its own at
##'   that instant. The receiver may shed it in turn, and pass it on down its
##'   own flows. The weights of the flows out of one pool sum to at most 1;
##'   what they leave over leaves the landscape.
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
flow <- function(from, to, type, out_of = 1L, into = out_of, weight = 1, ...) {
  kinds <- flow_kinds()
  if (missing(type) || !is.character(type) || length(type) != 1L ||
      !(type %in% names(kinds))) {
    stop("`type` is the kind of flow: one of ",
         paste0("\"", names(kinds), "\"", collapse = ", "), call. = FALSE)
  }
  ## A value no kind asked for is refused, so a misspelt argument cannot
  ## become a column nobody reads.
  extra <- list(...)
  own <- kinds[[type]]$columns
  if (length(extra) > 0L &&
      (is.null(names(extra)) || !all(names(extra) %in% own))) {
    stop("A ", type, " flow takes ",
         if (length(own) == 0L) "no further values" else paste(own, collapse = ", "),
         "; it was given ",
         paste(setdiff(names(extra) %||% "an unnamed value", own), collapse = ", "),
         call. = FALSE)
  }
  if (length(from) == 0L && length(to) == 0L) {
    return(empty_flows())
  }
  n <- max(length(from), length(to), length(out_of), length(into),
           length(weight), lengths(extra))
  recycle <- function(x, what) {
    if (!(length(x) %in% c(1L, n))) {
      stop("'", what, "' must be length 1 or ", n, call. = FALSE)
    }
    rep_len(x, n)
  }
  if (!is_index(out_of) || !is_index(into)) {
    stop("`out_of` and `into` are indices, counted from 1", call. = FALSE)
  }
  flows <- new_flows(recycle(from, "from"), recycle(to, "to"), type,
                     recycle(as.integer(out_of), "out_of"),
                     recycle(as.integer(into), "into"),
                     recycle(weight, "weight"))
  for (column in names(extra)) {
    flows[[column]] <- recycle(extra[[column]], column)
  }
  check_flows(flows)
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
    arriving <- passed[passed$to == patch, , drop = FALSE]
    arriving <- arriving[order(arriving$time), , drop = FALSE]
    ## An error here is the patch's, so it is reported under its name.
    results[[patch]] <- tryCatch({
      sc <- x$patches[[patch]]
      for (kind in unique(arriving$type)) {
        sc <- kinds[[kind]]$arriving(sc, arriving[arriving$type == kind, , drop = FALSE])
      }
      run(sc, record = record, refine_schedule = refine_schedule,
          method = method)
    }, error = function(e) {
      stop("Patch \"", patch, "\"",
           if (nrow(arriving) > 0L) ", with what flowed into it",
           ": ", conditionMessage(e), call. = FALSE)
    })
    passed <- rbind(passed, transfers_out_of(patch, results[[patch]], x$flows, kinds))
  }

  ## In time order; transfers at one instant keep the order their patches ran in.
  passed <- passed[order(passed$time), , drop = FALSE]
  rownames(passed) <- NULL
  structure(list(landscape = x, patches = results[names(x$patches)],
                 transfers = passed, method = method),
            class = "plant_landscape_run")
}

## Every amount passed between patches: when, from where to where, and how
## much reached the receiver.
empty_transfers <- function() {
  data.frame(time = numeric(0), from = character(0), to = character(0),
             type = character(0), out_of = integer(0), into = integer(0),
             amount = numeric(0), stringsAsFactors = FALSE)
}

## What a patch that has run passes down each of its flows: the flow's weight
## times each amount its kind says left by it.
transfers_out_of <- function(patch, result, flows, kinds) {
  flows <- flows[flows$from == patch, , drop = FALSE]
  out <- lapply(seq_len(nrow(flows)), function(i) {
    flow <- flows[i, , drop = FALSE]
    left <- tryCatch(kinds[[flow$type]]$leaving(result, flow), error = function(e) {
      stop("Patch \"", patch, "\", in what left it by a ", flow$type,
           " flow: ", conditionMessage(e), call. = FALSE)
    })
    if (!is.data.frame(left) || !all(c("time", "amount") %in% names(left)) ||
        !all(is.finite(left$time)) || !all(is.finite(left$amount)) ||
        any(left$amount < 0)) {
      stop("Patch \"", patch, "\" reported what left it by a ", flow$type,
           " flow as something other than finite times and finite, ",
           "non-negative amounts", call. = FALSE)
    }
    n <- nrow(left)
    data.frame(time = left$time, from = rep(patch, n), to = rep(flow$to, n),
               type = rep(flow$type, n), out_of = rep(flow$out_of, n),
               into = rep(flow$into, n), amount = flow$weight * left$amount,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, c(list(empty_transfers()), out))
}

##' @export
print.plant_landscape_run <- function(x, ...) {
  cat(sprintf("<plant_landscape_run> %d patch(es)%s, %d transfer(s) between them\n",
              length(x$patches),
              if (identical(x$method, "stochastic")) " (stochastic)" else "",
              nrow(x$transfers)))
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
