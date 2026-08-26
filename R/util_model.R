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
  ## Refuse a name that is not a parameter of this model (#636).
  ##
  ## ⚠️ WITHOUT THIS, A MISSPELT TRAIT IS A SILENT NO-OP AND THE RUN LOOKS FINE.
  ## The assignment below is `pars[trait_names] <- xi` on a list, which APPENDS an
  ## unknown name rather than dropping it -- so the junk element reaches the
  ## generated `Rcpp::as<*_Pars>`, which reads only the fields it knows and never
  ## looks at the extras (`RcppR6_post.hpp`: "No current support for a hook").
  ## The strategy then builds, solves and reports at the DEFAULT value of the
  ## parameter the caller believed they were varying. A trait sweep written with a
  ## typo runs the whole sweep at one point; a calibration fits a model it never
  ## perturbed. There is nothing in the output to say so, which is why this has to
  ## be refused rather than warned about.
  ##
  ## ⚠️ AFTER hyperpar(), NOT BEFORE. The hyperpar is what turns input traits into
  ## the parameter columns actually assigned -- it consumes e.g. `K_s` and emits
  ## `stem_P50`/`stem_c`, and its own guard already refuses an input that collides
  ## with something it derives. Checking the caller's original column names would
  ## reject every legitimate hyperparameterised trait.
  ##
  ## ⚠️ THE VALID SET IS `pars` ONLY, and that is narrower than "settable on the
  ## strategy". TF24f's `k_acclim`, `psi_fd_step` and `use_ad_gradient` are real
  ## user-settable fields living at the TOP level, and the assignment below cannot
  ## reach them, so they are not traits today. This turns that silence into an
  ## error, which is the honest report; making them reachable is a feature and is
  ## deliberately not done here (#636).
  known <- names(strategy$pars)
  ## Traits the hyperpar reads from the matrix but never assigns to the strategy.
  ## It declares them itself (see `make_FF16_hyperpar`), because `pars` cannot
  ## show them: FF16's `narea` derives `a_p1`, `a_p2` and `r_l` and is not a
  ## parameter, so a `narea` sweep is legitimate and must not be refused. A
  ## caller-supplied hyperpar with no declaration gets `pars` alone, which is the
  ## honest default -- we cannot know what an arbitrary closure reads.
  consumed <- attr(hyperpar, "input_traits", exact = TRUE)
  if (is.null(consumed)) consumed <- character(0)

  unknown <- setdiff(trait_names, c(known, consumed))
  if (length(unknown) > 0L) {
    ## Name the nearest valid parameter for each offender: the whole failure mode
    ## is a name that is almost right, so "did you mean" is the message doing the
    ## real work.
    candidates <- c(known, consumed)
    suggest <- vapply(unknown, function(u) {
      d <- utils::adist(u, candidates, ignore.case = TRUE)[1, ]
      near <- candidates[d == min(d) & d <= max(2L, ceiling(nchar(u) / 3))]
      if (length(near) == 0L) u else sprintf("%s (did you mean %s?)", u,
                                             paste(near, collapse = " / "))
    }, character(1L), USE.NAMES = FALSE)
    ## The strategy carries no `name` field at R level, so the class is what names
    ## the model -- e.g. "TF24_Strategy". `strategy$name` silently yields NULL.
    model <- class(strategy)[[1L]]
    stop(sprintf(paste0(
      "Unknown trait name%s for %s: %s\n",
      "A trait must name a parameter in the strategy's `pars`, or a trait the ",
      "hyperpar consumes; see names(%s()$pars) for the %d available.\n",
      "Note the hyperpar's own arguments (B_kl1, B_ks1, ...) are NOT traits: ",
      "pass them to make_%s_hyperpar() instead."),
      if (length(unknown) > 1L) "s" else "",
      model,
      paste(suggest, collapse = ", "),
      model, length(known), sub("_Strategy$", "", model)),
      call. = FALSE)
  }

  ## Assign only the columns that ARE parameters. A consumed trait like `narea`
  ## has already done its work inside the hyperpar, and appending it to `pars`
  ## would leave a junk element that the generated `Rcpp::as<*_Pars>` ignores --
  ## harmless to the model, but it makes a strategy compare unequal to a default
  ## one for a reason that has nothing to do with the model.
  ##
  ## ⚠️ `xi` below is UNNAMED and positional, aligned with `trait_names`, so this
  ## has to carry the positions rather than the names. `xi[assign_names]` returns
  ## a vector of NAs -- silently, and every parameter would land as NA.
  assign_at <- which(trait_names %in% known)
  assign_names <- trait_names[assign_at]

  f <- function(xi, br) {
    # Every column produced by hyperpar() is a biological parameter, which now
    # lives in the nested `pars` sub-object. Copy-back form: active/nested list
    # access returns a copy, so modify it then assign it back.
    pars <- strategy$pars
    pars[assign_names] <- xi[assign_at]
    strategy$pars <- pars
    if (is.list(br)) {
      if (!is.numeric(br$x) || length(br$x) < 2) {
        stop("birth_rate$x must be a numeric vector of at least length 2")
      }
      if (br$x[[1]] != 0) {
        stop(sprintf(
          "birth_rate$x must start at 0, not %g.\n",
          br$x[[1]]
        ))
      }
      max_lifetime <- p$max_patch_lifetime
      if (!is.null(max_lifetime) && max(br$x) < max_lifetime) {
        stop(sprintf(paste0(
          "birth_rate$x must extend to at least max_patch_lifetime (%g), ",
          "but ends at %g.\n",
          "Extend birth_rate$x (and birth_rate$y) to cover the full simulation duration."
        ), max_lifetime, max(br$x)))
      }
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
