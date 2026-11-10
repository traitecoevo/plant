##' Convert a matrix of trait values into a list of \code{Strategy} objects,
##' one per row, by applying the hyperparameter function and inserting the
##' resulting parameters into a copy of the default strategy.
##'
##' @title Generate strategies from traits
##' @param p A \code{Parameters} object containing a default strategy to
##'   modify. The trade-off function it carries in \code{p$hyperpar} (see
##'   \code{\link{default_parameters}}) is applied unless \code{hyperpar} is
##'   given.
##' @param traits Trait values as a \emph{matrix}, with column names
##'   corresponding to traits (see \code{\link{trait_matrix}}); one row per
##'   strategy. A column may name a field of the strategy's \code{pars}, a
##'   trait the trade-off function reads, or a numeric top-level strategy
##'   field (TF24t's thermal parameters); \code{\link{parameter_table}} lists
##'   them.
##' @param hyperpar Trade-off (hyperparameter) function to use: a closure
##'   \code{f(m, s, filter = TRUE)} that translates ecological traits (e.g.
##'   \code{lma}, wood density) into the low-level strategy parameters,
##'   encoding the model's trade-offs. Defaults to the one \code{p} carries
##'   (\code{param_hyperpar(p)}); \code{NULL} is the explicit identity, which
##'   writes every trait straight to the strategy. A \code{p} carrying none is
##'   treated as the identity too, and said so once per session.
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
  hyperpar <- resolve_hyperpar(hyperpar)
  traits <- hyperpar(traits, strategy)

  trait_names <- colnames(traits)
  known <- names(strategy$pars)
  fields <- strategy_field_names(strategy)
  refuse_unknown_traits(trait_names, strategy, known, fields, hyperpar)

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
  field_at <- which(trait_names %in% fields)
  field_names <- trait_names[field_at]

  f <- function(xi, br) {
    # Every column produced by hyperpar() that is a biological parameter lives
    # in the nested `pars` sub-object. Copy-back form: active/nested list
    # access returns a copy, so modify it then assign it back.
    pars <- strategy$pars
    pars[assign_names] <- xi[assign_at]
    strategy$pars <- pars
    # A top-level strategy field (TF24t's thermal parameters) is assigned where
    # it lives, beside `pars`.
    for (j in seq_along(field_at)) {
      strategy[[field_names[[j]]]] <- xi[[field_at[[j]]]]
    }
    set_strategy_birth_rate(strategy, br, p$max_patch_lifetime)
  }

  # insert custom traits and birth values into default strategy template
  mapply(f, matrix_to_list(traits), birth_rate, SIMPLIFY = FALSE)
}

# The trade-off function generate_strategy() applies: `NULL` is the explicit
# identity, and a function is used as given.
resolve_hyperpar <- function(hyperpar) {
  ## `hyperpar = NULL` is the explicit identity. A `p` with no trade-offs
  ## attached gives the same identity through param_hyperpar(), but flagged, so
  ## that a caller who expected a trade-off (TF24's K_s -> stem_P50, say) hears
  ## that none ran. Once per session: this is information, not an error, and a
  ## sweep would otherwise say it a thousand times.
  if (is.null(hyperpar)) {
    hyperpar <- identity_hyperpar
  }
  if (!is.function(hyperpar)) {
    stop("`hyperpar` must be a trade-off function f(m, s, filter = TRUE), or NULL",
         call. = FALSE)
  }
  if (isTRUE(attr(hyperpar, "unattached", exact = TRUE))) {
    rlang::inform(paste0(
      "No trade-off function is attached to these Parameters, so traits are ",
      "written straight to the strategy. default_parameters(model) attaches ",
      "plant's testing trade-offs; the fieldtrial package has the recommended ",
      "ones and curated species/site sets. This message is shown once per session."),
      .frequency = "once", .frequency_id = "plant_no_tradeoffs_attached")
  }
  hyperpar
}

# Stops unless every trait column names something a strategy of this model
# takes: a field of `pars`, a numeric top-level field, or a trait the trade-off
# function declares it consumes.
refuse_unknown_traits <- function(trait_names, strategy, known, fields, hyperpar) {
  ## Refuse a name that is not a parameter of this model (#636).
  ##
  ## ⚠️ WITHOUT THIS, A MISSPELT TRAIT IS A SILENT NO-OP AND THE RUN LOOKS FINE.
  ## generate_strategy() assigns traits into `pars` as a list, which APPENDS an
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
  ## The valid set is `pars` plus the numeric top-level strategy fields (TF24t's
  ## thermal parameters, TF24f's `k_acclim`), each assigned where it lives by
  ## generate_strategy(). `strategy_field_names()` keeps only numeric scalars, so a field of
  ## any other type stays unreachable from a trait matrix and is refused here
  ## rather than coerced.
  ## Traits the hyperpar reads from the matrix but never assigns to the strategy.
  ## It declares them itself (see `frozen_FF16_hyperpar`), because `pars` cannot
  ## show them: FF16's `narea` derives `a_p1`, `a_p2` and `r_l` and is not a
  ## parameter, so a `narea` sweep is legitimate and must not be refused. A
  ## caller-supplied hyperpar with no declaration gets `pars` alone, which is the
  ## honest default -- we cannot know what an arbitrary closure reads.
  consumed <- attr(hyperpar, "input_traits", exact = TRUE)
  if (is.null(consumed)) consumed <- character(0)

  unknown <- setdiff(trait_names, c(known, fields, consumed))
  if (length(unknown) > 0L) {
    ## Name the nearest valid parameter for each offender: the whole failure mode
    ## is a name that is almost right, so "did you mean" is the message doing the
    ## real work.
    candidates <- c(known, fields, consumed)
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
      "A trait must name a parameter in the strategy's `pars`, a numeric ",
      "top-level strategy field, or a trait the trade-off function consumes; ",
      "parameter_table(\"%s\") lists the %d available.\n",
      "Note the trade-off function's own coefficients (B_kl1, B_ks1, ...) are ",
      "NOT traits: they are set where the trade-off is built."),
      if (length(unknown) > 1L) "s" else "",
      model,
      paste(suggest, collapse = ", "),
      sub("_Strategy$", "", model), length(c(known, fields, consumed))),
      call. = FALSE)
  }
  invisible(NULL)
}

# `strategy` with the birth rate `br` set: a constant, or a list of control
# points `x` (from 0 to at least `max_lifetime`) and values `y`.
set_strategy_birth_rate <- function(strategy, br, max_lifetime) {
  if (is.list(br)) {
    if (!is.numeric(br$x) || length(br$x) < 2) {
      stop("birth_rate$x must be a numeric vector of at least length 2", call. = FALSE)
    }
    if (br$x[[1]] != 0) {
      stop(sprintf(
        "birth_rate$x must start at 0, not %g.\n",
        br$x[[1]]
      ), call. = FALSE)
    }
    if (!is.null(max_lifetime) && max(br$x) < max_lifetime) {
      stop(sprintf(paste0(
        "birth_rate$x must extend to at least max_patch_lifetime (%g), ",
        "but ends at %g.\n",
        "Extend birth_rate$x (and birth_rate$y) to cover the full simulation duration."
      ), max_lifetime, max(br$x)), call. = FALSE)
    }
    strategy$birth_rate_x <- br$x
    strategy$birth_rate_y <- br$y
    strategy$is_variable_birth_rate <- TRUE
  } else if (is.numeric(br)) {
    strategy$birth_rate_y <- br
    strategy$is_variable_birth_rate <- FALSE
  } else {
    stop("Invalid type in birth_rate - need either a list with x, y control points or a numeric", call. = FALSE)
  }
  strategy
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
##' @param hyperpar Trade-off (hyperparameter) function to use. Defaults to the
##'   one \code{p} carries (\code{param_hyperpar(p)}, see
##'   \code{\link{default_parameters}}); \code{NULL} is the explicit identity.
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

