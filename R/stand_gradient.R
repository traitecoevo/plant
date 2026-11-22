# The model's census entry points, each bound to the model's name as the
# deterministic runner's own class gives it. The C++ side dispatches on that name
# and refuses a model that declares no census; TF24t is refused here, by name,
# because the reason is not that no census was written.
census_backend <- function(solver) {
  model <- extract_RcppR6_template_types(solver, "DeterministicPatchRunner")[[1]]
  if (identical(model, "TF24t")) {
    stop("No census is declared for the TF24t model, so it has no stand ",
         "gradient: plant does not put the thermal-damage layer's states on ",
         "the tape, and a gradient without them would describe TF24 rather ",
         "than TF24t")
  }
  list(metric_names = function(x) census_metric_names(x, model),
       census = function(x) census_values(x, model),
       state_adjoint = function(x) census_state_adjoint(x, model),
       trait_names = function(x) census_trait_names(x, model),
       undifferentiable = function() census_undifferentiable(model),
       trait_gradient = function(x, metrics) census_trait_gradient(x, model, metrics),
       trait_gradient_forward = function(x, columns, metrics)
         census_trait_gradient_forward(x, model, columns, metrics),
       control = function(x) census_gradient_control(x, model))
}

# The models the forward route walks. TF24 and its subclasses solve the leaf
# inside every rate evaluation and keep that solve's choices, which a forward
# walk would solve again rather than replay.
forward_gradient_models <- "FF16"

# Below this many trait columns, mode = "auto" takes the forward route. A forward
# walk costs about a run plus a share per column, and does not depend on the
# metrics; the sweep costs the same whatever the number of columns, plus a
# share per metric. Measured on FF16 at lifetime 30, for an invader's single
# fitness metric, the forward route is the cheaper up to 8 columns, one walk at
# 8 lanes, and a ninth column needs a second walk, which costs more than the
# sweep. A resident, or more metrics, moves the break-even further out, so the
# invader's single metric sets it.
forward_gradient_crossover <- 9L

# Whether to take the forward route, for `n` trait columns of `solver`.
use_forward_gradient <- function(solver, mode, n) {
  switch(mode,
         reverse = FALSE,
         forward = TRUE,
         auto = is.null(forward_gradient_refusal(solver)) &&
           n < forward_gradient_crossover)
}

# Why the forward route cannot take `solver`'s run, or NULL where it can. The
# C++ entry refuses the same runs; this says so before anything runs.
forward_gradient_refusal <- function(solver) {
  model <- extract_RcppR6_template_types(solver, "DeterministicPatchRunner")[[1]]
  if (!model %in% forward_gradient_models) {
    return(paste0(
      "it is available for ", paste(forward_gradient_models, collapse = ", "),
      " only, not the ", model, " model",
      if (startsWith(model, "TF24")) {
        paste0(": its strategy keeps the choices of its leaf solve, which a ",
               "forward walk would solve again rather than replay")
      }))
  }
  NULL
}

# mode = "forward" on a run the forward route cannot take is refused before
# anything runs, saying why.
check_gradient_mode <- function(solver, mode) {
  if (identical(mode, "forward")) {
    why <- forward_gradient_refusal(solver)
    if (!is.null(why)) {
      # The sweep is the alternative only for a model that declares a census.
      model <- extract_RcppR6_template_types(solver, "DeterministicPatchRunner")[[1]]
      stop("mode = \"forward\" refused: ", why,
           if (model %in% c("FF16", "TF24", "TF24f")) ". Use mode = \"reverse\"")
    }
  }
}

# Every metric the run reports, named, with each one's kind ("census" or
# "fitness") beside it. The C++ side lists them in one order for the value, the
# tangent and the sweep, and this is the one place R reads that list.
stand_metric_table <- function(solver) {
  b <- census_backend(solver)
  named <- b$metric_names(solver)
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
##' differentiable (FF16, TF24 and TF24f), each species' net reproduction ratio and offspring
##' production, named \code{"<species>.net_reproduction_ratio"} and
##' \code{"<species>.offspring_production"}. They are the same numbers as
##' a run's \code{net_reproduction_ratios} and \code{offspring_production},
##' reported under the names \code{\link{stand_gradient}} differentiates them
##' by.
##'
##' @param x A run \code{DeterministicPatchRunner} object, or a \code{plant_run}, for the FF16, TF24 or TF24f model, on the
##'   birth-date coordinate (TF24's default; \code{control(node_density_coordinate = "birth_date")} for FF16).
##' @return A named numeric vector, one entry per metric.
##' @export
stand_census <- function(x) {
  x <- as_solver(x)
  stand_metrics_of_kind(x, "census")
}

##' @rdname stand_census
##' @export
stand_fitness <- function(x) {
  x <- as_solver(x)
  stand_metrics_of_kind(x, "fitness")
}

stand_metrics_of_kind <- function(solver, kind) {
  b <- census_backend(solver)
  table <- stand_metric_table(solver)
  # An empty vector reads as a stand with no species, and indexing it by name
  # fails far from the cause; say instead that the model declares none.
  if (!any(table$kind == kind)) {
    stop("This model declares no differentiable ", kind, " metric")
  }
  values <- stats::setNames(b$census(solver), table$name)
  values[table$kind == kind]
}

##' Sensitivity of the stand metrics to the ODE state at the current time.
##'
##' One row per metric -- the census metrics, then for a model that declares
##' its fitness differentiable each species' net reproduction ratio and offspring
##' production, named as \code{\link{stand_gradient}} names them -- and one
##' column per ODE state entry, in the order
##' \code{x$patch$ode_state} writes. This is what the reverse pass starts
##' with. On the birth-date coordinate, the only one the gradient runs on, the
##' trapezium weights are constants, so the cohort-height columns carry the
##' integrand's response alone.
##'
##' @param x A run \code{DeterministicPatchRunner} object, or a \code{plant_run}, for the FF16, TF24 or TF24f model, on the
##'   birth-date coordinate (TF24's default; \code{control(node_density_coordinate = "birth_date")} for FF16).
##' @return A numeric matrix, metrics by ODE state entries.
##' @export
stand_census_state_adjoint <- function(x) {
  x <- as_solver(x)
  b <- census_backend(x)
  rows <- b$state_adjoint(x)
  out <- do.call(rbind, rows)
  rownames(out) <- stand_metric_table(x)$name
  out
}

##' The \code{Control} entries a census gradient depends on.
##'
##' The entries inside a rate evaluation and the schedule that move the
##' gradient. Three are tolerances of a model's inner solves: the leaf's
##' (\code{GSS_tol_abs}, \code{ci_abs_tol}) and a newborn's height
##' (\code{offspring_height_tol}). One is a finite-difference step
##' (\code{node_gradient_eps}) and one the schedule's refinement tolerance
##' (\code{schedule_eps}). The sixth, \code{gradient_curvature_floor}, moves no
##' forward number and still decides which rows exist, by refusing a collar
##' response the profit curvature is too small to support. The list is the same
##' for every model, so it holds entries a model does not read (FF16 reads
##' neither leaf tolerance nor the floor): two of its gradients that differ only
##' there are refused although they are equal, and no two that differ in one
##' are compared. The ODE integrator's tolerances and step bounds move the
##' trajectory too, and are the run's own \code{Control}, not listed here.
##' \code{stand_gradient} records these and \code{stand_gradient_compare}
##' refuses a pair that disagrees.
##'
##' @param x A run \code{DeterministicPatchRunner} object, or a \code{plant_run}, for the FF16, TF24 or TF24f model, on the
##'   birth-date coordinate (TF24's default; \code{control(node_density_coordinate = "birth_date")} for FF16).
##' @return A named numeric vector of length six.
##' @export
gradient_control <- function(x) {
  x <- as_solver(x)
  census_backend(x)$control(x)
}

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
##' @param x A run \code{DeterministicPatchRunner} object, or a \code{plant_run}, for the FF16, TF24 or TF24f model, on the
##'   birth-date coordinate (TF24's default; \code{control(node_density_coordinate = "birth_date")} for FF16).
##'   The sweep reads the state at every step. A run that kept none (one run
##'   without \code{record = "states"}) is run again inside this call to keep
##'   them, which costs as much as the run did and gives the same numbers; the
##'   solver keeps them afterwards (\code{x$keep_states} is then \code{TRUE}),
##'   so a second call does not repeat it. Run with \code{record = "states"} to
##'   pay for the run once.
##' @param metrics Metric names to differentiate; defaults to every metric the
##'   model declares: the census metrics (\code{\link{stand_census}}) and, for
##'   a model that declares its fitness differentiable, each species'
##'   \code{"<species>.net_reproduction_ratio"} and
##'   \code{"<species>.offspring_production"} (\code{\link{stand_fitness}}).
##'   The steps are recorded once for all the metrics asked for, and each
##'   further metric adds only a walk back over them, a small part of the
##'   first metric's cost; under the sweep, asking for fewer traits saves
##'   nothing (see \code{mode}).
##' @param traits Column names to differentiate with respect to; defaults to
##'   every differentiable trait every species declares. A column is named for
##'   its species and its parameter, as \code{"1.lma"}.
##' @param hyperpar The hyperparameter function the strategies were built with,
##'   by default the one \code{add_strategies()} uses; or \code{NULL} for
##'   partial derivatives. See the section on traits and parameters.
##' @param mode How the gradient is computed. \code{"reverse"} sweeps back over
##'   the recorded steps once per metric, at the same cost whatever the number
##'   of traits. \code{"forward"} walks forward over the same steps carrying one
##'   direction per trait column, up to eight in a walk, so its cost grows with
##'   the columns and does not depend on the number of metrics. It is for FF16
##'   only, and an error for any other model. \code{"auto"},
##'   the default, takes the forward route where it can when fewer than nine
##'   columns are needed, and the sweep otherwise. The two routes give the same
##'   shape, names and refusals; a partial derivative agrees to rounding (about
##'   1e-15 of the metric's largest entry), and a total through a
##'   hyperparameter function can lose more to cancellation in the chain rule.
##'   Under a hyperparameter function the columns counted are the ones the
##'   chain rule reads: the traits asked for and each parameter derived from
##'   them that moves with one (FF16's lma adds \code{k_l} and \code{r_l}).
##' @return A list with \code{value}, the metrics at the end of the run;
##'   \code{gradient}, a metrics-by-traits matrix; \code{refusal}, one entry per
##'   metric holding the reason and the species it was found in (\code{NA}
##'   where it belongs to none), or \code{NULL} where the metric answered; and
##'   \code{control}, the entries the gradient was taken at.
##'
##'   Refused outright, with an error: a run on the height coordinate or in a
##'   constant field; with the absorbing top layer on, a run with a harvest or
##'   climate extreme that reached any cohort; a strategy no gradient exists
##'   for whatever its run (TF24's deep-crown shading or sapwood acclimation,
##'   \code{a_sw > 0}), refused before a run that kept no states is repeated;
##'   and a trade-off that derives a parameter an equation reads and that has
##'   no column.
##'
##'   A run resumed from \code{\link{export_patch_state}()} is differentiated
##'   from its resume time, with the exported state held fixed: how that state
##'   depends on the traits is not in the gradient. For an invader of a resumed
##'   resident that is exact, since the invader's own cohorts all start after
##'   the resume. For the resident it is an approximation, to be stated beside
##'   any resident gradient taken that way. With the absorbing top layer on
##'   (\code{Control$absorb_epsilon}), which cohorts are absorbed is held as the
##'   run chose it, and what they carry into the layer and the height it is
##'   pinned at move with the traits.
##'
##' @section Traits and parameters:
##' A column is the derivative as its trait moves through the construction that
##' built the strategy, which is what a difference through \code{add_strategies()}
##' measures. The hyperparameter function derives some parameters from the traits
##' (TF24's testing trade-offs, \code{default_hyperpar("TF24")}, which sets twelve parameters, among them \code{k_l}, \code{r_l}, \code{d_I},
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
##' With \code{hyperpar = NULL} every parameter the model can differentiate
##' has a column, and each is a partial derivative with every other parameter
##' fixed.
##'
##' A column can be zero because of another parameter's value. FF16, TF24
##' and TF24f offer the coarse-root pool's \code{a_cr1}, \code{r_cr} and
##' \code{k_cr}. At the default \code{a_cr1 = 0} the pool is empty, so the
##' columns of its respiration and turnover rates, \code{r_cr} and
##' \code{k_cr}, are exact zeros by construction, not findings about the pool;
##' \code{a_cr1}'s own column is live there. Set \code{a_cr1 > 0} to read
##' the other two.
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
##' After \code{x$run_mutant(p)}, a gradient is of the invaders' metrics in
##' the residents' recorded field, which the invaders read at each rate
##' evaluation and do not move, so it carries no derivative: the rare-mutant
##' derivative. A run of \code{x$run_mutant_constant()} is refused.
##' An invasion that kept no states is run again, as an invasion, to keep them;
##' set \code{x$keep_states <- TRUE} before it to save that run.
##' Several invaders in \code{p} share one run and each has its own columns,
##' and the sweep's cost grows with their cohorts, so a set of invaders costs
##' about what running them one at a time does. An invader of a resumed resident starts at
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
##' Do not average such a grid, and do not form a covariance from one. A
##' single outlying cell dominates any second moment, which is the first thing an
##' active-subspace or sensitivity workflow computes, so the reduction inherits
##' the step placement instead of the ecology. Treat a gradient far outside its
##' neighbourhood as a diagnostic: nudge the trait a fraction of a percent and
##' take the answer that is stable, or resolve the event by forcing a step at the
##' crossing.
##' @export
stand_gradient <- function(x, metrics = NULL, traits = NULL,
                           hyperpar = param_hyperpar(x$parameters),
                           mode = c("auto", "reverse", "forward")) {
  x <- as_solver(x)
  mode <- match.arg(mode)
  check_gradient_mode(x, mode)
  ## The identity (hyperpar = NULL, or Parameters carrying no trade-offs)
  ## derives nothing, so the total is the partial.
  if (!is.null(hyperpar) && !is_identity_hyperpar(hyperpar)) {
    return(stand_gradient_through(x, metrics, traits, hyperpar, mode))
  }
  b <- census_backend(x)
  all_metrics <- stand_metric_table(x)$name
  # Every parameter the strategy carries has a column, bar the few the model
  # states it cannot differentiate. So the width does not depend on what the
  # sweep can currently answer, and a caller indexing by position sees the same
  # shape from run to run.
  all_traits <- b$trait_names(x)
  if (is.null(metrics)) {
    metrics <- all_metrics
  }
  if (length(metrics) == 0L) {
    stop("`metrics` names no metric; leave it NULL for every metric")
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
    # saying so sends readers looking for a spelling mistake. The model states why
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

  value <- stats::setNames(b$census(x), all_metrics)[metrics]
  # Only the metrics asked for are swept. Named, not positioned: C++ resolves
  # the names against the same list it reports them from.
  if (use_forward_gradient(x, mode, length(unique(traits)))) {
    # The forward route computes only the columns asked for, one lane each;
    # the sweep computes every column at once.
    columns <- unique(traits)
    swept <- b$trait_gradient_forward(x, match(columns, all_traits),
                                      as.character(metrics))
  } else {
    columns <- all_traits
    swept <- b$trait_gradient(x, as.character(metrics))
  }
  gradient <- matrix(unlist(swept$gradient), nrow = length(metrics),
                     ncol = length(columns), byrow = TRUE,
                     dimnames = list(metrics, columns))

  list(value = value,
       gradient = gradient[metrics, traits, drop = FALSE],
       refusal = stats::setNames(swept$refusal, metrics),
       control = gradient_control(x))
}

# The chain rule through a hyperparameter function: a trait's total is its partial
# plus G[, derived] %*% d(derived)/d(trait), with the derived partials from the
# sweep or the forward walk.
stand_gradient_through <- function(solver, metrics, traits, hyperpar,
                                   mode = "auto") {
  if (!is.function(hyperpar)) {
    stop("`hyperpar` must be a trade-off function f(m, s, filter = TRUE), as ",
         "param_hyperpar(p) or default_hyperpar(\"TF24\") give, or NULL for ",
         "partial derivatives")
  }
  strategies <- solver$parameters$strategies
  species <- seq_along(strategies)
  if (!is.null(traits)) {
    species <- species[species %in% sub("^([0-9]+)\\..*", "\\1", traits)]
  }
  trade_offs <- Map(trade_off, strategies[species], species,
                    MoreArgs = list(hyperpar = hyperpar))
  derived <- unlist(lapply(trade_offs, rownames))
  columns <- census_backend(solver)$trait_names(solver)
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
  unread <- attr(census_backend(solver)$undifferentiable(), "unread")
  lost <- setdiff(derived, read)
  lost <- lost[!trait_without_species(lost) %in% unread]
  if (length(lost) > 0L) {
    stop("`hyperpar` derives ", paste(lost, collapse = ", "), " from the ",
         "species' traits, which has no gradient column, so the chain rule ",
         "cannot be completed")
  }
  # A derived parameter that no trait asked for moves adds only terms of exactly
  # zero to the sum. The sweep computes every column at once and is given them
  # all; on the forward route each column is a lane, so it is given only the
  # derived parameters that move.
  moved <- unlist(lapply(trade_offs, function(J) {
    along <- intersect(colnames(J), traits)
    rownames(J)[rowSums(J[, along, drop = FALSE] != 0) > 0]
  }))
  wanted <- union(traits, intersect(read, moved))
  if (use_forward_gradient(solver, mode, length(wanted))) {
    mode <- "forward"
  } else {
    mode <- "reverse"
    wanted <- union(traits, read)
  }
  partial <- stand_gradient(solver, metrics, wanted, hyperpar = NULL,
                            mode = mode)
  G <- partial$gradient
  total <- G[, traits, drop = FALSE]
  for (J in trade_offs) {
    along <- intersect(colnames(J), traits)
    from <- intersect(rownames(J), colnames(G))
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
##' \code{Control}: five of the six entries change the trajectory and so change
##' the function being differentiated, and the sixth,
##' \code{gradient_curvature_floor}, decides which rows exist. See
##' \code{\link{gradient_control}}.
##'
##' @param a,b Results of \code{stand_gradient}.
##' @return The element-wise difference \code{a$gradient - b$gradient}, over the
##'   metrics and traits both carry.
##' @export
stand_gradient_compare <- function(a, b) {
  if (!setequal(names(a$control), names(b$control))) {
    stop("These gradients record different Control entries, so they were taken ",
         "on different models or builds: ",
         paste(union(setdiff(names(a$control), names(b$control)),
                     setdiff(names(b$control), names(a$control))),
               collapse = ", "))
  }
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
