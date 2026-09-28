##' Census metrics of a solved stand.
##'
##' The trapezium integral of \eqn{n_k \psi(\mathrm{state}_k)} over the size
##' distribution, for each metric \eqn{\psi}. The grid is the coordinate the
##' density is carried in -- the birth date, or the height -- and the inflow
##' boundary node closes it, so the interval between it and the smallest cohort
##' is in the sum.
##'
##' @param scm A run \code{SCM} object for the TF24 strategy.
##' @return A named numeric vector, one entry per census metric.
##' @export
stand_census <- function(scm) {
  stats::setNames(census_tf24(scm), census_metric_names_tf24())
}

##' Sensitivity of the census metrics to the ODE state at the current time.
##'
##' One row per metric and one column per ODE state entry, in the order
##' \code{scm$patch$ode_state} writes. This is what the reverse pass is seeded
##' with. The cohort-height columns carry the trapezium weights as well as the
##' integrand only where the height is the coordinate integrated over; on the
##' birth-date coordinate the weights are constants.
##'
##' @param scm A run \code{SCM} object for the TF24 strategy.
##' @return A numeric matrix, metrics by ODE state entries.
##' @export
stand_census_state_adjoint <- function(scm) {
  rows <- census_state_adjoint_tf24(scm)
  out <- do.call(rbind, rows)
  rownames(out) <- census_metric_names_tf24()
  out
}

##' The \code{Control} entries a census gradient depends on.
##'
##' The first four move the trajectory the gradient is taken along, so two
##' gradients taken at different values of any of them are gradients of
##' different functions. \code{gradient_curvature_floor} is here for a second
##' reason: it moves no forward number at all and still decides which rows
##' exist, by refusing a collar response the profit curvature is too small to
##' support. \code{stand_gradient} records them and
##' \code{stand_gradient_compare} refuses a pair that disagrees.
##'
##' @param scm A run \code{SCM} object for the TF24 strategy.
##' @return A named numeric vector of length five.
##' @export
gradient_control <- function(scm) gradient_control_tf24(scm)

# The parameter a gradient column names, with its species index stripped. An
# unknown name is matched on this so a missing species prefix is distinguishable
# from a parameter the strategy does not declare.
trait_without_species <- function(x) sub("^[0-9]+\\.", "", x)

##' Gradient of a stand's census metrics with respect to traits.
##'
##' Doubles in and doubles out: the active scalar is created and destroyed inside
##' one call and never crosses into R.
##'
##' Two terms, because a census reads the traits as well as the state. The reverse
##' pass is seeded on the states the census reads at the end of the run and runs
##' back over the steps the adaptive pass resolved; added to it is the census's own
##' reading of the traits at that state, which is not a sensitivity of the state
##' and so no sweep produces it.
##'
##' @param scm A run \code{SCM} object for the TF24 strategy.
##' @param metrics Census metric names to differentiate; defaults to all of them.
##' @param traits Column names to differentiate with respect to; defaults to
##'   every differentiable parameter every species declares. A column is named
##'   for its species and its parameter, as \code{"1.lma"}. With \code{hyperpar},
##'   these name the traits to differentiate through it, and must be given.
##' @param hyperpar \code{NULL}, or the hyperparameter function the stand was
##'   built with, as \code{TF24_hyperpar}. See the section on traits and
##'   parameters.
##' @return A list with \code{value}, the metrics at the end of the run;
##'   \code{gradient}, a metrics-by-traits matrix; \code{refusal}, one entry per
##'   metric holding the reason and where it was found, or \code{NULL} where the
##'   metric answered; and \code{control}, the entries the gradient was taken at.
##'   With \code{hyperpar}, also \code{jacobian}: for each trait, the derivative of
##'   every parameter the hyperparameter function derives from it.
##'
##' @section Traits and parameters:
##' ⚠️ **Without \code{hyperpar}, each column is a PARTIAL derivative with every
##' other parameter held fixed**, including those a hyperparameter function
##' derives from the trait. \code{TF24_hyperpar} sets \code{k_l}, \code{r_l} and
##' \code{nmass_l} from \code{lma}, so \code{gradient[, "1.lma"]} is not the
##' derivative a trait fit wants: on a 3-year stand it was 2.4 to 4 times a
##' finite difference of \code{lma} through \code{add_strategies()}.
##'
##' Pass \code{hyperpar} for the total derivative with respect to the traits as a
##' trait matrix gives them. Each trait's column is its own partial plus every
##' derived parameter's partial times that parameter's derivative with respect
##' to the trait. The hyperparameter function is closed-form, so that derivative
##' is taken by central difference of the function itself, costing no model run.
##' Pass the function the stand was built with: a run does not record it.
##'
##'   A refused metric's whole row is \code{NaN}: a sum has no defined value with
##'   an undefined term, so refusal is metric-level and carries no localisation
##'   within a metric. Every other number is one the sweep computed, an exact
##'   zero included — a parameter no gradient exists for cannot be asked for, so
##'   it never reaches a column. \code{psi_crit} and \code{root_psi_crit} are
##'   two such: each is derived from its curve's \code{(P50, c)} pair, which
##'   \code{phylloptim} derives for itself, so setting one reaches no equation
##'   and asking for its gradient errors rather than returning a zero.
##'
##' @section A single entry can dominate a second moment:
##' A gradient taken where the run crosses an unresolved event — a cohort
##' reaching \code{hmat} between two steps, say — is a derivative of the
##' trajectory's step placement rather than of the model, and it can be two
##' orders larger than its neighbours. One measured cell was 211 times the
##' median of its grid.
##'
##' ⚠️ **Do not average such a grid, and do not form a covariance from one.** A
##' single outlying cell dominates any second moment, which is the first thing an
##' active-subspace or sensitivity workflow computes, so the reduction inherits
##' the step placement instead of the ecology. Treat a gradient far outside its
##' neighbourhood as a diagnostic: nudge the trait a fraction of a percent and
##' take the answer that is stable, or resolve the event by forcing a step at the
##' crossing.
##' @export
stand_gradient <- function(scm, metrics = NULL, traits = NULL, hyperpar = NULL) {
  if (!is.null(hyperpar)) {
    return(stand_gradient_through(scm, metrics, traits, hyperpar))
  }
  all_metrics <- census_metric_names_tf24()
  # Every parameter the strategy carries has a column, bar the few the model
  # states it cannot differentiate. So the width does not depend on what the
  # sweep can currently answer, and a caller indexing by position sees the same
  # shape from run to run.
  all_traits <- census_trait_names_tf24(scm)
  if (is.null(metrics)) {
    metrics <- all_metrics
  }
  if (is.null(traits)) {
    traits <- all_traits
  }
  unknown <- setdiff(metrics, all_metrics)
  if (length(unknown) > 0L) {
    stop("Unknown census metric: ", paste(unknown, collapse = ", "))
  }
  unknown <- setdiff(traits, all_traits)
  if (length(unknown) > 0L) {
    # A parameter the model carries but cannot differentiate is not unknown, and
    # saying so sent readers looking for a spelling mistake. The model states why
    # it has no gradient; that sentence is the refusal.
    why <- census_undifferentiable_tf24()
    named <- trait_without_species(unknown)
    cannot <- named %in% names(why)
    if (any(cannot)) {
      stop("No gradient exists for ",
           paste(sprintf("%s (%s)", unknown[cannot], why[named[cannot]]),
                 collapse = "; "),
           ". Ask for other traits.")
    }
    # A bare parameter name is the shape of the defect the prefix exists to
    # prevent, so say what the columns are called rather than only refusing.
    hint <- if (any(named %in% trait_without_species(all_traits))) {
      paste0(". Columns carry their species index, as \"",
             all_traits[[1]], "\"")
    } else {
      ""
    }
    stop("Unknown trait: ", paste(unknown, collapse = ", "), hint)
  }

  value <- stand_census(scm)[metrics]
  # Only the metrics asked for are swept. A metric costs a sweep of the whole
  # trajectory, so computing all three and subsetting the answer charged a
  # caller who wanted one for three. Named, not positioned: C++ resolves the
  # names against the same list it reports them from.
  swept <- census_trait_gradient_tf24(scm, as.character(metrics))
  gradient <- do.call(rbind, swept$gradient)
  dimnames(gradient) <- list(metrics, all_traits)

  list(value = value,
       gradient = gradient[metrics, traits, drop = FALSE],
       refusal = stats::setNames(swept$refusal, metrics),
       control = gradient_control(scm))
}

# The chain rule through a hyperparameter function. A trait the caller fits (lma,
# rho, K_s...) reaches the model both as itself and through every parameter the
# hyperpar derives from it (lma sets k_l, r_l and nmass_l), and stand_gradient's
# columns are partials with the derived ones held fixed. The total is
#     d(metric)/d(trait) = G[, trait] + sum_k G[, k] * d(p_k)/d(trait)
# with d(p_k)/d(trait) from a central difference of the hyperpar itself. The
# hyperpar is closed-form algebra, so that difference is accurate to ~1e-10 and
# costs no model evaluation.
stand_gradient_through <- function(scm, metrics, traits, hyperpar) {
  if (!is.function(hyperpar)) {
    stop("`hyperpar` must be a hyperparameter function, as TF24_hyperpar")
  }
  if (is.null(traits) || length(traits) == 0L) {
    stop("Name the traits to differentiate through `hyperpar`, as \"1.lma\"")
  }
  partial <- stand_gradient(scm, metrics)
  G <- partial$gradient
  strategies <- scm$parameters$strategies
  undifferentiable <- names(census_undifferentiable_tf24())

  total <- matrix(NA_real_, nrow(G), length(traits),
                  dimnames = list(rownames(G), traits))
  jacobian <- list()
  for (tr in traits) {
    species <- suppressWarnings(as.integer(sub("\\..*$", "", tr)))
    name <- trait_without_species(tr)
    if (is.na(species) || species < 1L || species > length(strategies) ||
        identical(name, tr)) {
      stop("A trait names its species and its column, as \"1.lma\": got \"",
           tr, "\"")
    }
    s <- strategies[[species]]
    x <- s$pars[[name]]
    if (is.null(x)) {
      stop("\"", tr, "\" is not a parameter of species ", species,
           "'s strategy, so its value cannot be read to differentiate at")
    }
    derived <- function(v) {
      out <- hyperpar(trait_matrix(v, name), s, filter = FALSE)
      out[1, setdiff(colnames(out), name), drop = TRUE]
    }
    h <- 6e-6 * max(abs(x), 1e-3)
    dp <- (derived(x + h) - derived(x - h)) / (2 * h)
    dp <- dp[dp != 0]
    cols <- paste0(species, ".", names(dp))
    # A derived parameter with no column contributes nothing only if no equation
    # reads it, which the model states; anything else leaves the sum incomplete.
    missing <- !(cols %in% colnames(G))
    unread <- names(dp) %in% undifferentiable
    if (any(missing & !unread)) {
      stop("`hyperpar` derives ", paste(names(dp)[missing & !unread],
                                        collapse = ", "),
           " from ", name, ", which has no gradient column, so the chain rule ",
           "for \"", tr, "\" cannot be completed")
    }
    keep <- !missing
    own <- if (tr %in% colnames(G)) G[, tr] else 0
    total[, tr] <- own + G[, cols[keep], drop = FALSE] %*% dp[keep]
    jacobian[[tr]] <- dp
  }

  list(value = partial$value,
       gradient = total,
       refusal = partial$refusal,
       control = partial$control,
       jacobian = jacobian)
}

##' Which of a gradient's metrics were refused.
##'
##' Refusal is metric-level: a sum has no defined value with an undefined term,
##' so a refused metric's whole row is \code{NaN} and one reason serves it.
##'
##' @param g A result of \code{stand_gradient}.
##' @return A named logical, one entry per metric.
##' @export
stand_gradient_refused <- function(g) {
  vapply(g$refusal, function(r) !is.null(r), logical(1))
}

##' Compare two stand gradients.
##'
##' Two gradients are comparable only if they were taken at the same
##' \code{Control}: four of the five entries change the trajectory and so change
##' the function being differentiated, and the fifth,
##' \code{gradient_curvature_floor}, decides which rows exist. See
##' \code{\link{gradient_control}}.
##'
##' @param a,b Results of \code{stand_gradient}.
##' @return The element-wise difference \code{a$gradient - b$gradient}, over the
##'   metrics and traits both carry.
##' @export
stand_gradient_compare <- function(a, b) {
  differing <- names(a$control)[!identical_doubles(a$control, b$control)]
  if (length(differing) > 0L) {
    stop("These gradients were taken at different Control values, so they are ",
         "gradients of different functions: ",
         paste(differing, collapse = ", "))
  }
  metrics <- intersect(rownames(a$gradient), rownames(b$gradient))
  traits <- intersect(colnames(a$gradient), colnames(b$gradient))
  a$gradient[metrics, traits, drop = FALSE] -
    b$gradient[metrics, traits, drop = FALSE]
}

# Element-wise exact equality, NA-safe, for two named numeric vectors of the
# same names.
identical_doubles <- function(x, y) {
  vapply(names(x), function(n) identical(unname(x[[n]]), unname(y[[n]])),
         logical(1))
}
