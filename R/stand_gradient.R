# The model's census entry points, chosen off the SCM's own class. Every
# exported function below reads them from here, so a model that declares no
# census is refused once, by name, instead of failing on a C++ type check.
census_backend <- function(scm) {
  type <- extract_RcppR6_template_types(scm, "SCM")[[1]]
  switch(type,
    TF24 = list(metric_names = census_metric_names_tf24,
                census = census_tf24,
                state_adjoint = census_state_adjoint_tf24,
                trait_names = census_trait_names_tf24,
                undifferentiable = census_undifferentiable_tf24,
                trait_gradient = census_trait_gradient_tf24,
                control = gradient_control_tf24),
    FF16 = list(metric_names = census_metric_names_ff16,
                census = census_ff16,
                state_adjoint = census_state_adjoint_ff16,
                trait_names = census_trait_names_ff16,
                undifferentiable = census_undifferentiable_ff16,
                trait_gradient = census_trait_gradient_ff16,
                control = gradient_control_ff16),
    # Named rather than left to the default, because the reason is not that no
    # census was written: no derivative of the thermal-damage layer exists.
    TF24t = stop("No census is declared for the TF24t model, so it has no ",
                 "stand gradient: phylloptim refuses the derivative of its ",
                 "leaf with the thermal-damage layer on, and a gradient ",
                 "without it would describe TF24 rather than TF24t"),
    stop("No census is declared for the ", type, " model, so it has no ",
         "stand gradient; TF24 and FF16 declare one"))
}

# Every metric the run reports, named, with each one's kind ("census" or
# "fitness") beside it. The C++ side lists them in one order for the value, the
# tangent and the sweep, and this is the one place R reads that list.
stand_metric_table <- function(scm) {
  b <- census_backend(scm)
  named <- b$metric_names(scm)
  list(name = unname(named), kind = names(named))
}

##' Metrics of a solved stand.
##'
##' \code{stand_census} reports the census metrics: the trapezium integral of
##' \eqn{n_k \psi(\mathrm{state}_k)} over the size distribution, for each
##' per-plant quantity \eqn{\psi}, summed over the species. The grid is the
##' coordinate the density is carried in -- the birth date, or the height -- and
##' the inflow boundary node closes it, so the interval between it and the
##' smallest cohort is in the sum.
##'
##' \code{stand_fitness} reports, for a model that declares its fitness
##' differentiable (FF16 and TF24), each species' net reproduction ratio and offspring
##' production, named \code{"<species>.net_reproduction_ratio"} and
##' \code{"<species>.offspring_production"}. They are the same numbers as
##' \code{scm$net_reproduction_ratios} and \code{scm$offspring_production},
##' reported under the names \code{\link{stand_gradient}} differentiates them
##' by.
##'
##' @param scm A run \code{SCM} object for the TF24 or FF16 model, on the
##'   birth-date coordinate (TF24's default; \code{control(node_density_coordinate = "birth_date")} for FF16).
##' @return A named numeric vector, one entry per metric.
##' @export
stand_census <- function(scm) {
  stand_metrics_of_kind(scm, "census")
}

##' @rdname stand_census
##' @export
stand_fitness <- function(scm) {
  stand_metrics_of_kind(scm, "fitness")
}

stand_metrics_of_kind <- function(scm, kind) {
  b <- census_backend(scm)
  table <- stand_metric_table(scm)
  # An empty vector reads as a stand with no species, and indexing it by name
  # fails far from the cause; say instead that the model declares none.
  if (!any(table$kind == kind)) {
    stop("This model declares no differentiable ", kind, " metric")
  }
  values <- stats::setNames(b$census(scm), table$name)
  values[table$kind == kind]
}

##' Sensitivity of the stand metrics to the ODE state at the current time.
##'
##' One row per metric -- the census metrics, then for a model that declares
##' its fitness differentiable each species' net reproduction ratio and offspring
##' production, named as \code{\link{stand_gradient}} names them -- and one
##' column per ODE state entry, in the order
##' \code{scm$patch$ode_state} writes. This is what the reverse pass starts
##' with. The cohort-height columns carry the trapezium weights as well as the
##' integrand only where the height is the coordinate integrated over; on the
##' birth-date coordinate the weights are constants.
##'
##' @param scm A run \code{SCM} object for the TF24 or FF16 model, on the
##'   birth-date coordinate (TF24's default; \code{control(node_density_coordinate = "birth_date")} for FF16).
##' @return A numeric matrix, metrics by ODE state entries.
##' @export
stand_census_state_adjoint <- function(scm) {
  b <- census_backend(scm)
  rows <- b$state_adjoint(scm)
  out <- do.call(rbind, rows)
  rownames(out) <- stand_metric_table(scm)$name
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
##' @param scm A run \code{SCM} object for the TF24 or FF16 model, on the
##'   birth-date coordinate (TF24's default; \code{control(node_density_coordinate = "birth_date")} for FF16).
##' @return A named numeric vector of length five.
##' @export
gradient_control <- function(scm) census_backend(scm)$control(scm)

# The parameter a gradient column names, with its species index stripped. An
# unknown name is matched on this so a missing species prefix is distinguishable
# from a parameter the strategy does not declare.
trait_without_species <- function(x) sub("^[0-9]+\\.", "", x)

##' Gradient of a stand's metrics with respect to traits.
##'
##' Doubles in and doubles out: the active scalar is created and destroyed inside
##' one call and never crosses into R.
##'
##' Two terms, because a census reads the traits as well as the state. The reverse
##' pass starts on the states the census reads at the end of the run and runs
##' back over the steps the adaptive pass resolved; added to it is the census's own
##' reading of the traits at that state, which is not a sensitivity of the state
##' and so no sweep produces it.
##'
##' @param scm A run \code{SCM} object for the TF24 or FF16 model, on the
##'   birth-date coordinate (TF24's default; \code{control(node_density_coordinate = "birth_date")} for FF16).
##' @param metrics Metric names to differentiate; defaults to every metric the
##'   model declares: the census metrics (\code{\link{stand_census}}) and, for
##'   a model that declares its fitness differentiable, each species'
##'   \code{"<species>.net_reproduction_ratio"} and
##'   \code{"<species>.offspring_production"} (\code{\link{stand_fitness}}).
##'   Each metric costs one sweep of the recorded trajectory, so name the ones
##'   wanted.
##' @param traits Column names to differentiate with respect to; defaults to
##'   every differentiable trait every species declares. A column is named for
##'   its species and its parameter, as \code{"1.lma"}.
##' @param hyperpar The hyperparameter function the strategies were built with,
##'   by default the one \code{add_strategies()} uses; or \code{NULL} for
##'   partial derivatives. See the section on traits and parameters.
##' @return A list with \code{value}, the metrics at the end of the run;
##'   \code{gradient}, a metrics-by-traits matrix; \code{refusal}, one entry per
##'   metric holding the reason and where it was found, or \code{NULL} where the
##'   metric answered; and \code{control}, the entries the gradient was taken at.
##'
##' @section Traits and parameters:
##' A column is the derivative as its trait moves through the construction that
##' built the strategy, which is what a difference through \code{add_strategies()}
##' measures. The hyperparameter function derives some parameters from the traits
##' (\code{TF24_hyperpar} sets 13, among them \code{k_l}, \code{r_l}, \code{d_I},
##' \code{a_f3} and the stem curve, from six traits including \code{lma}), so a
##' trait's column is its own partial plus each derived parameter's partial times
##' that parameter's derivative with respect to the trait. The function is
##' closed-form, so that derivative is a central difference of the function
##' itself, and costs no model run.
##'
##' Every derived parameter leaves the columns, as
##' \code{colnames(hyperpar(matrix(numeric(0), 1, 0), s, filter = FALSE))} lists
##' them for a strategy \code{s}. Each requested species' derived parameters are
##' checked against its traits first: a strategy built by hand (TF24's default is
##' not on its trade-off), edited after it was built, built at an input it does not
##' keep (FF16's \code{narea}), or built through another function, is refused
##' rather than differentiated along a trade-off it does not lie on. Pass the
##' function it was built with, or \code{hyperpar = NULL}.
##'
##' With \code{hyperpar = NULL} every parameter has a column, and each is a
##' partial derivative with every other parameter fixed.
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
##' @section An invader against a recorded resident:
##' After \code{scm$run_mutant(p)}, a gradient is of the invaders' metrics in
##' the residents' recorded field, which the invaders read at each rate
##' evaluation and do not move, so it carries no derivative: the rare-mutant
##' derivative. A run of \code{scm$run_mutant_constant()} is refused.
##' An invasion that kept no states is run again, as an invasion, to keep them;
##' set \code{scm$keep_states <- TRUE} before it to save that run.
##' Several invaders in \code{p} share one run and each has its own columns,
##' but each metric is still one sweep, so a set of invaders costs about what
##' running them one at a time does. An invader of a resumed resident starts at
##' the resume time, as \code{run_mutant} requires.
##'
##' @section The birth rate:
##' For a species run at a constant birth rate the column
##' \code{"<species>.birth_rate"} is the derivative with respect to that rate,
##' which reaches the density every cohort is born at and the offspring
##' integral. A species with a variable birth rate has no such column.
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
stand_gradient <- function(scm, metrics = NULL, traits = NULL,
                           hyperpar = param_hyperpar(scm$parameters)) {
  if (!is.null(hyperpar)) {
    return(stand_gradient_through(scm, metrics, traits, hyperpar))
  }
  b <- census_backend(scm)
  all_metrics <- stand_metric_table(scm)$name
  # Every parameter the strategy carries has a column, bar the few the model
  # states it cannot differentiate. So the width does not depend on what the
  # sweep can currently answer, and a caller indexing by position sees the same
  # shape from run to run.
  all_traits <- b$trait_names(scm)
  if (is.null(metrics)) {
    metrics <- all_metrics
  }
  if (is.null(traits)) {
    traits <- all_traits
  }
  unknown <- setdiff(metrics, all_metrics)
  if (length(unknown) > 0L) {
    stop("Unknown metric: ", paste(unknown, collapse = ", "))
  }
  unknown <- setdiff(traits, all_traits)
  if (length(unknown) > 0L) {
    # A parameter the model carries but cannot differentiate is not unknown, and
    # saying so sent readers looking for a spelling mistake. The model states why
    # it has no gradient; that sentence is the refusal.
    why <- b$undifferentiable()
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

  value <- stats::setNames(b$census(scm), all_metrics)[metrics]
  # Only the metrics asked for are swept. A metric costs a sweep of the whole
  # trajectory, so computing all three and subsetting the answer charged a
  # caller who wanted one for three. Named, not positioned: C++ resolves the
  # names against the same list it reports them from.
  swept <- b$trait_gradient(scm, as.character(metrics))
  gradient <- do.call(rbind, swept$gradient)
  dimnames(gradient) <- list(metrics, all_traits)

  list(value = value,
       gradient = gradient[metrics, traits, drop = FALSE],
       refusal = stats::setNames(swept$refusal, metrics),
       control = gradient_control(scm))
}

# The chain rule through a hyperparameter function: a trait's total is its partial
# plus G[, derived] %*% d(derived)/d(trait), with the derived partials from the sweep.
stand_gradient_through <- function(scm, metrics, traits, hyperpar) {
  if (!is.function(hyperpar)) {
    stop("`hyperpar` must be a hyperparameter function, as TF24_hyperpar, or ",
         "NULL for partial derivatives")
  }
  strategies <- scm$parameters$strategies
  species <- seq_along(strategies)
  if (!is.null(traits)) {
    species <- species[species %in% sub("^([0-9]+)\\..*", "\\1", traits)]
  }
  trade_offs <- Map(trade_off, strategies[species], species,
                    MoreArgs = list(hyperpar = hyperpar))
  derived <- unlist(lapply(trade_offs, rownames))
  columns <- census_backend(scm)$trait_names(scm)
  if (is.null(traits)) {
    traits <- setdiff(columns, derived)
  }
  named <- intersect(traits, derived)
  if (length(named) > 0L) {
    stop(paste(named, collapse = ", "), " derived from the species' traits by ",
         "`hyperpar`: ask for the trait it comes from, or pass hyperpar = NULL ",
         "for its partial derivative")
  }
  # A derived parameter without a column adds nothing only where the model states
  # no equation reads it; any other leaves the sum incomplete.
  read <- intersect(derived, columns)
  unread <- names(census_backend(scm)$undifferentiable())
  lost <- setdiff(derived, read)
  lost <- lost[!trait_without_species(lost) %in% unread]
  if (length(lost) > 0L) {
    stop("`hyperpar` derives ", paste(lost, collapse = ", "), " from the ",
         "species' traits, which has no gradient column, so the chain rule ",
         "cannot be completed")
  }
  partial <- stand_gradient(scm, metrics, union(traits, read), hyperpar = NULL)
  G <- partial$gradient
  total <- G[, traits, drop = FALSE]
  for (J in trade_offs) {
    along <- intersect(colnames(J), traits)
    from <- intersect(rownames(J), read)
    total[, along] <- total[, along] +
      G[, from, drop = FALSE] %*% J[from, along, drop = FALSE]
  }
  partial$gradient <- total
  partial
}

# Each parameter `hyperpar` derives for species `i`, differentiated in each trait
# by one call on the traits and each moved up and down: rows and columns named.
trade_off <- function(s, i, hyperpar) {
  derived <- colnames(hyperpar(matrix(numeric(0), 1L, 0L), s, filter = FALSE))
  if (length(derived) == 0L) {
    return(matrix(numeric(0), 0L, 0L))
  }
  x <- unlist(s$pars[setdiff(names(s$pars), derived)])
  P <- length(x)
  h <- 6e-6 * pmax(abs(x), 1e-3)
  at <- matrix(x, P, P, byrow = TRUE)
  rows <- rbind(x, at + diag(h, P), at - diag(h, P))
  colnames(rows) <- names(x)
  out <- hyperpar(rows, s, filter = FALSE)[, derived, drop = FALSE]
  # ⚠️ CHECKED: a chain through a trade-off the strategy is not on returns finite,
  # wrong numbers. 1e-7, as the function's filter keeps a default to 1.5e-8.
  gap <- abs(out[1L, ] - unlist(s$pars)[derived])
  off <- !is.finite(gap) | gap > 1e-7 * abs(out[1L, ])
  if (any(off)) {
    unkept <- attr(hyperpar, "input_traits", exact = TRUE)
    stop("Species ", i, "'s ", paste(derived[off], collapse = ", "),
         " are not what `hyperpar` gives at its traits: it was built by hand, ",
         "edited after it was built",
         if (length(unkept)) paste0(", built at a ", paste(unkept, collapse = ", "),
                                    " the strategy does not keep"),
         ", or built through another function. Pass the function it was built ",
         "with, or hyperpar = NULL for partial derivatives")
  }
  up <- out[1L + seq_len(P), , drop = FALSE]
  down <- out[1L + P + seq_len(P), , drop = FALSE]
  J <- t(up - down) / rep(2 * h, each = length(derived))
  dimnames(J) <- list(paste0(i, ".", derived), paste0(i, ".", names(x)))
  J[, colSums(J != 0) > 0, drop = FALSE]
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
