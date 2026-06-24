##' Convert a matrix of trait values into a list of \code{Strategy} objects,
##' one per row, by applying the hyperparameter function and inserting the
##' resulting parameters into a copy of the default strategy.
##'
##' @title Generate strategies from traits
##' @param p A \code{Parameters} object containing a default strategy to
##'   modify.  Any hyperparameterisation included will be applied.
##' @param traits Trait values as a \emph{matrix}, with column names
##'   corresponding to traits (see \code{\link{trait_matrix}}); one row per
##'   strategy.
##' @param hyperpar Hyperparameter function to use. By default links to the
##'   standard function for this strategy type. It translates ecological traits
##'   (e.g. \code{lma}, wood density) into the low-level strategy parameters,
##'   encoding the model's trade-offs.
##' @param birth_rate Birth rate(s) for each row of \code{traits}: a scalar or
##'   vector (constant birth rate, set as \code{strategy$birth_rate_y}) or a
##'   list with \code{x}, \code{y} control points (a varying birth rate, which
##'   also sets \code{strategy$birth_rate_x} and
##'   \code{is_variable_birth_rate = TRUE}).
##' @param x Deprecated (\code{strategy_list}); use \code{traits}.
##' @param parameters Deprecated (\code{strategy_list}); use \code{p}.
##' @param birth_rate_list Deprecated; use \code{birth_rate}.
##'
##' @export
##' @rdname generate_strategy
generate_strategy <- function(p, traits, hyperpar = param_hyperpar(p),
                              birth_rate = 1) {
  if (!is.matrix(traits)) {
    stop("Invalid type traits -- expected a matrix")
  }

  strategy <- p$strategy_default
  traits <- hyperpar(traits, strategy)

  trait_names <- colnames(traits)
  f <- function(xi, br) {
    # Every column produced by hyperpar() is a biological parameter, which now
    # lives in the nested `pars` sub-object. Copy-back form: active/nested list
    # access returns a copy, so modify it then assign it back.
    pars <- strategy$pars
    pars[trait_names] <- xi
    strategy$pars <- pars
    if (is.list(br)) {
      strategy$birth_rate_x <- br$x
      strategy$birth_rate_y <- br$y
      strategy$is_variable_birth_rate <- TRUE
    } else if (is.numeric(br)) {
      strategy$birth_rate_y <- br
      strategy$is_variable_birth_rate <- FALSE
    } else {
      stop("Invalid type in birth_rate - need either a list with x, y control points or a numeric")
    }
    strategy
  }

  # insert custom traits and birth values into default strategy template
  mapply(f, matrix_to_list(traits), birth_rate, SIMPLIFY = FALSE)
}

##' Helper function to create trait matrices suitable for
##' \code{\link{generate_strategy}} and \code{\link{add_strategies}}.
##'
##' @title Create trait matrix
##' @param x Values
##' @param trait_name Name of a single trait
##' @export
##' @author Rich FitzJohn
trait_matrix <- function(x, trait_name) {
  m <- matrix(x, ncol=length(trait_name))
  colnames(m) <- trait_name
  m
}

##' Add strategies to a \code{Parameters} object. \code{add_strategies} appends
##' the new strategies to any existing residents (controlled by
##' \code{keep_existing}); \code{add_mutant} introduces strategies as mutants
##' (i.e. replacing the resident set). Both translate trait values into
##' strategies via \code{\link{generate_strategy}} and are pipe-friendly
##' (the \code{Parameters} object is the first argument).
##'
##' @title Add strategies (or a mutant) to Parameters
##' @param p A \code{Parameters} object.
##' @param traits A matrix of traits corresponding to the new strategies to
##'   introduce (see \code{\link{trait_matrix}}).
##' @param hyperpar Hyperparameter function to use. By default links to the
##'   standard function for this strategy type.
##' @param birth_rate Birth rate(s), one per row of \code{traits}. See
##'   \code{\link{generate_strategy}}.
##' @param keep_existing Should existing resident strategies be retained?
##' @param trait_matrix Deprecated (\code{expand_parameters}/\code{mutant_parameters});
##'   use \code{traits}.
##' @param birth_rate_list Deprecated; use \code{birth_rate}.
##' @param keep_existing_strategies Deprecated; use \code{keep_existing}.
##' @export
##' @rdname add_strategies
add_strategies <- function(p, traits, hyperpar = param_hyperpar(p),
                           birth_rate = 1, keep_existing = TRUE) {

  if (nrow(traits) != length(birth_rate)) {
    stop("Must provide exactly one birth rate input for each species")
  }
  extra <- generate_strategy(p, traits, hyperpar, birth_rate)
  n_extra <- length(extra)

  ret <- p <- validate(p) # Ensure times are set up correctly.

  ## Determine node introduction times
  if (length(p$strategies) == 0L) {
    times_new <- p$node_schedule_times_default
  } else {
    ## if residents are present, use all unique times of all residents
    times_new <- unique(sort(unlist(p$node_schedule_times)))
  }

  if (keep_existing) {
    ret$strategies <- c(p$strategies, extra)
    ret$node_schedule_times <- c(p$node_schedule_times,
                                  rep(list(times_new), n_extra))
  } else {
    ret$strategies <- extra
    ret$node_schedule_times <- rep(list(times_new), n_extra)
  }

  ## Clear this if it's present:
  attr(ret, "offspring_production") <- NULL

  ret
}

##' @export
##' @rdname add_strategies
add_mutant <- function(p, traits, hyperpar = param_hyperpar(p),
                       birth_rate = 1) {
  add_strategies(p, traits, hyperpar = hyperpar, birth_rate = birth_rate,
                 keep_existing = FALSE)
}

##' @export
##' @rdname generate_strategy
strategy_list <- function(x, parameters, hyperpar = param_hyperpar(parameters),
                          birth_rate_list = 1) {
  .Deprecated("generate_strategy")
  generate_strategy(parameters, x, hyperpar = hyperpar,
                    birth_rate = birth_rate_list)
}

##' @export
##' @rdname add_strategies
expand_parameters <- function(trait_matrix, p, hyperpar = param_hyperpar(p),
                              birth_rate_list = 1,
                              keep_existing_strategies = TRUE) {
  .Deprecated("add_strategies")
  add_strategies(p, trait_matrix, hyperpar = hyperpar,
                 birth_rate = birth_rate_list,
                 keep_existing = keep_existing_strategies)
}

##' @export
##' @rdname add_strategies
mutant_parameters <- function(trait_matrix, p, hyperpar = param_hyperpar(p),
                              birth_rate_list = 1,
                              keep_existing_strategies = FALSE) {
  .Deprecated("add_mutant")
  add_strategies(p, trait_matrix, hyperpar = hyperpar,
                 birth_rate = birth_rate_list,
                 keep_existing = keep_existing_strategies)
}

remove_residents <- function(p) {
  if (length(p$strategies) > 0L) {
    p$strategies <- list()
    p$birth_rate <- numeric(0)
    p$node_schedule_times <- list()
  }
  p
}
