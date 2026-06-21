##' Construct a \code{Control} object. \code{control()} is a lowercase alias for
##' the \code{Control()} constructor, whose defaults are the pragmatic,
##' fast-ish settings used for essentially all of plant's runs (see
##' \code{control.cpp}). \code{control_accurate()} tightens the ODE and schedule
##' tolerances for high-accuracy runs at the cost of speed.
##'
##' @title Control presets
##' @param ... Named control fields, passed to \code{Control()}.
##' @param base An optional \code{Control} object to tighten; defaults are used
##'   if omitted.
##' @return A \code{Control} object.
##' @rdname control_presets
##' @export
control <- function(...) Control(...)

##' @rdname control_presets
##' @export
control_accurate <- function(base = Control()) {
  base$ode_tol_rel       <- 1e-6
  base$ode_tol_abs       <- 1e-6
  base$ode_step_size_max <- 1e-1
  base$schedule_eps      <- 1e-3
  base
}


##' Basic default settings for a given strategy, environment only really
##' used for templating initially and will be overloaded later by passing
##' an environment to the SCM API (suggesting perhaps the template could be
##' removed).
##' @title Basic default parameters for a given strategy
##' @author Rich FitzJohn
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @param env And environment object
##' @export
scm_base_parameters <- function(type = NA, env = environment_type(type)) {
  Parameters(type, env)()
}


##' Run the SCM.
##'
##' The node-introduction schedule can be adaptively refined in C++ by setting
##' \code{refine_schedule = TRUE} (this replaces the former \code{build_schedule}
##' function). Setting \code{collect = TRUE} returns tidied output collected at
##' every ODE step (replacing the former \code{run_scm_collect}); otherwise the
##' \code{SCM} object itself is returned for interrogation.
##'
##' @title Run SCM
##' @param p Parameters object
##' @param env Environment object (defaults to the strategy's environment)
##' @param ctrl Control object
##' @param refine_schedule Should the node-introduction schedule be adaptively
##'   refined before/while running (using \code{schedule_eps} and
##'   \code{schedule_nsteps} from \code{ctrl})?
##' @param collect Should tidied results be collected at every step and
##'   returned (instead of the \code{SCM} object)?
##' @param use_ode_times Should ODE times be used?
##' @return When \code{collect = FALSE}, an \code{SCM} object. When
##'   \code{collect = TRUE}, a list of tidied patch output with
##'   \code{offspring_production}, \code{net_reproduction_ratios} and the
##'   (possibly refined) parameters \code{p}.
##' @author Rich FitzJohn
##' @rdname run_scm
##' @export
run_scm <- function(p, env = NULL,
                    ctrl = control(),
                    refine_schedule = FALSE, collect = FALSE,
                    use_ode_times = FALSE) {

  types <- extract_RcppR6_template_types(p, "Parameters")

  if (is.null(env))
    env <- Environment(types[[1]])

  scm <- do.call('SCM', types)(p, env, ctrl)
  if (use_ode_times) {
    # Pin integration to the schedule's ode_times (loaded from p$ode_times).
    sched <- scm$node_schedule
    sched$use_ode_times <- TRUE
    scm$node_schedule <- sched
  }
  if (collect) {
    scm$collect <- TRUE
  }

  if (refine_schedule) {
    scm$refine_schedule()
  } else {
    scm$run()
  }

  if (!collect) {
    return(scm)
  }

  results <- lapply(scm$history, "[[", "state") |> tidy_patch()
  results[["offspring_production"]] <- scm$offspring_production
  results[["net_reproduction_ratios"]] <- scm$net_reproduction_ratios
  results[["p"]] <- scm$parameters

  results
}

