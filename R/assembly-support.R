##' Helper function for creating parameter objects suitable for an
##' assembly.
##' @title Helper function for creating parameter objects
##' @param ... Named set of parameters
##' @param pars A list of parameters
##' @param base_parameters_fn Function for creating base parameter set (default scm_base_parameters)
##' @param base_control_fn Function for creating base Control object (default scm_base_control)
##' @param make_hyperpar_fn Function for creating hyperparameterisation (default make_FF16_hyperpar)
##' @export
assembly_parameters <- function(..., pars = NULL, type = NA,
                                base_parameters_fn = scm_base_parameters,
                                base_control_fn = scm_base_control,
                                make_hyperpar_fn = make_FF16_hyperpar) {
  p <- base_parameters_fn(type)
  ctrl <- base_control_fn()

  ## These are nice to have:
  ctrl$equilibrium_solver_name <- "hybrid"
  ctrl$equilibrium_nsteps <- 60

  if (is.null(pars)) {
    pars <- list(...)
  } else if (length(list(...)) > 0L) {
    stop("Do not provide both ... and pars")
  }

  if (length(pars) > 0L) {
    assert_named_if_not_empty(pars)

    excl <- c("strategy_default", "hyperpar")
    pos <- setdiff(
      c(
        names(formals(make_hyperpar_fn)),
        names(p),
        names(p$strategy_default)
      ),
      excl
    )
    unk <- setdiff(names(pars), pos)
    if (length(unk) > 0L) {
      stop("Unknown parameters: ", paste(unk, collapse = ", "))
    }

    nms_hyper <- intersect(names(pars), names(formals(make_hyperpar_fn)))
    p <- modify_list(p, pars)
    p$strategy_default <- modify_list(p$strategy_default, pars)
  }
  p
}


##' Check low-abundance strategies for viability.
##'
##' @title Check low-abundance strategies for viability
##' @param p A Parameters object
##' @param ctrl Control object
##' @export
check_inviable <- function(p, ctrl) {
  ## eps_test: *Relative* value to use for determining what
  ## "low abundance" means.  Species that have a offspring arrival of less than
  ## `eps_test * max(p$birth_rate)` will be tested.  By default
  ##  this is 1 100th of the maximum offspring arrival.
  ## TODO: don't do anything if we don't have at least 2 species?
  eps <- ctrl$equilibrium_extinct_birth_rate
  ## TODO: This was ctrl$equilibrium_inviable_test, but I think
  ## that birth offspring arrival actually makes more sense?  It's fractional
  ## though so who knows.
  eps_test <- 1e-2
  birth_rate <- sapply(p$strategies, function(s) s$birth_rate_y, simplify = TRUE)
  ## NOTE: We don't actually run to equilibrium here; this is just
  ## because it's a useful way of doing incoming -> outgoing offspring
  ## rain.
  runner <- make_equilibrium_runner(p, ctrl =ctrl)
  offspring_production <- runner(birth_rate)
  
  test <- which(offspring_production < birth_rate &
                  birth_rate < max(offspring_production) * eps_test)
  test <- test[order(offspring_production[test])]
  
  drop <- logical(length(offspring_production))
  
  for (i in test) {
    plant_log_inviable(paste("Testing species", i),
                       stage="testing", species=i)
    x <- offspring_production
    x[i] <- eps
    res <- runner(x)
    if (res[[i]] < eps) {
      plant_log_inviable(paste("Removing species", i),
                         stage="removing", species=i)
      drop[[i]] <- TRUE
      res[[i]] <- 0.0
      offspring_production <- res
    }
  }
  
  ## It's possible that things slip through and get driven extinct by
  ## the time that they reach here.
  drop <- drop | offspring_production < eps
  
  attr(offspring_production, "drop") <- drop
  offspring_production
}
