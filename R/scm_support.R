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

##' Export the full state of a patch from a (run) \code{SCM} so it can be
##' re-imported to seed a new run (see \code{\link{set_initial_state}}). The
##' exported object captures everything needed to reproduce the patch's forward
##' trajectory: every node's ODE state plus the per-node birth bookkeeping
##' (introduction time, patch-age density and survival probability at birth)
##' that is not part of the ODE state but feeds the rates and lifetime-fitness
##' integrals, the patch age, and the not-yet-introduced ("residual") portion of
##' the node-introduction schedule.
##'
##' @title Export patch state from an SCM
##' @param scm An \code{SCM} object that has been run.
##' @param step Optional 1-based index into \code{scm$history} (requires the run
##'   to have been performed with \code{collect = TRUE}); when \code{NULL}
##'   (default) the SCM's current/final patch is used.
##' @return A list describing the patch state, suitable for
##'   \code{\link{set_initial_state}}: \code{time}, \code{n} (nodes per species),
##'   \code{ode_state} (flat), the per-species lists \code{node_times},
##'   \code{patch_density} and \code{pr_patch_survival}, and the residual
##'   \code{node_schedule_times}.
##' @seealso \code{\link{set_initial_state}}, \code{\link{run_scm}}
##' @export
export_patch_state <- function(scm, step = NULL) {
  patch <- if (is.null(step)) scm$patch else scm$history[[step]]
  time <- patch$time
  species <- patch$species

  ## Residual schedule: the introductions not yet represented among the seeded
  ## nodes. A patch snapshot sits integrated *to* its next introduction time
  ## without having introduced that node yet (its node times are all strictly
  ## below `time`), so the residual keeps every original time at or after `time`
  ## -- the resumed run introduces the node due at `time` itself, then the rest.
  tol <- 1e-8
  residual <- lapply(scm$parameters$node_schedule_times,
                     function(tt) tt[tt >= time - tol])

  list(
    time = time,
    n = as.integer(vapply(species, function(s) s$size, numeric(1))),
    ode_state = patch$ode_state,
    node_times = lapply(species, function(s) s$node_times),
    patch_density = lapply(species, function(s) s$patch_densities),
    pr_patch_survival = lapply(species, function(s) s$pr_patch_survival_at_birth),
    node_schedule_times = residual
  )
}

##' Build an initial patch state from a specified size distribution, for seeding
##' a patch at age 0 with pre-existing plants instead of growing it from empty
##' (the ecological motivation of \code{plant} issue #304). The returned object
##' is consumed by \code{\link{set_initial_state}}.
##'
##' Each species' initial nodes are described by their heights and densities;
##' the remaining ODE state (mortality, fecundity, accumulated reproduction,
##' heartwood, ...) starts at zero, all nodes are introduced at patch age 0, and
##' the recruitment schedule continues for \code{t > 0} (the seeded distribution
##' replaces the \code{t = 0} recruit). Pathologically large/dense initial
##' conditions can produce non-finite densities; \code{run_scm} guards against
##' this and errors with a suggestion to use more plausible inputs.
##'
##' @title Build an initial size distribution
##' @param p A \code{Parameters} object.
##' @param heights Per-species node heights: a numeric vector (single species) or
##'   a list of numeric vectors (one per species).
##' @param densities Per-species node densities (same shape as \code{heights}).
##'   Supply this or \code{log_densities}.
##' @param log_densities Per-species node log-densities (alternative to
##'   \code{densities}).
##' @param env Environment object (defaults to the strategy's environment).
##' @param ctrl Control object.
##' @return A state list suitable for \code{\link{set_initial_state}}.
##' @seealso \code{\link{set_initial_state}}, \code{\link{export_patch_state}}
##' @export
make_initial_state <- function(p, heights, densities = NULL,
                               log_densities = NULL, env = NULL,
                               ctrl = control()) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  if (is.null(env)) {
    env <- Environment(types[[1]])
  }
  n_spp <- length(p$strategies)

  as_list <- function(z) if (is.list(z)) z else list(z)
  heights <- as_list(heights)
  if (length(heights) != n_spp) {
    stop("`heights` must have one entry per species (", n_spp, ")")
  }
  if (is.null(log_densities)) {
    if (is.null(densities)) {
      stop("supply either `densities` or `log_densities`")
    }
    log_densities <- lapply(as_list(densities), log)
  } else {
    log_densities <- as_list(log_densities)
  }

  ## Learn the node ODE layout and the (node-free) environment ODE tail from a
  ## fresh patch, plus the patch-age-0 disturbance weights for birth bookkeeping.
  patch <- do.call("Patch", types)(p, env, ctrl)
  ode_names <- patch$species[[1]]$new_node$ode_names
  hi <- match("height", ode_names)
  ldi <- match("log_density", ode_names)
  if (is.na(hi) || is.na(ldi)) {
    stop("could not locate 'height'/'log_density' in node ODE names")
  }
  node_ode_size <- length(ode_names)
  env_state <- patch$ode_state # fresh patch has no nodes: environment ODE only
  pr_surv0 <- patch$pr_survival(0)
  dens0 <- patch$density(0)

  n <- integer(n_spp)
  ode_chunks <- vector("list", n_spp)
  node_times <- vector("list", n_spp)
  patch_density <- vector("list", n_spp)
  pr_patch_survival <- vector("list", n_spp)
  for (i in seq_len(n_spp)) {
    h <- heights[[i]]
    ld <- log_densities[[i]]
    if (length(h) != length(ld)) {
      stop("heights and densities must have equal length (species ", i, ")")
    }
    ## The model requires nodes ordered by decreasing height.
    o <- order(h, decreasing = TRUE)
    h <- h[o]
    ld <- ld[o]
    mat <- matrix(0, nrow = node_ode_size, ncol = length(h))
    mat[hi, ] <- h
    mat[ldi, ] <- ld
    ode_chunks[[i]] <- as.vector(mat) # column-major: node-by-node, matching set_ode_state
    n[i] <- length(h)
    node_times[[i]] <- rep(0, length(h))
    patch_density[[i]] <- rep(dens0, length(h))
    pr_patch_survival[[i]] <- rep(pr_surv0, length(h))
  }

  list(
    time = 0,
    n = n,
    ode_state = c(unlist(ode_chunks, use.names = FALSE), env_state),
    node_times = node_times,
    patch_density = patch_density,
    pr_patch_survival = pr_patch_survival,
    ## Continue recruitment for t > 0; the seeded distribution replaces the t=0 recruit.
    node_schedule_times = lapply(p$node_schedule_times, function(tt) tt[tt > 1e-8])
  )
}

##' Write an exported patch state (from \code{\link{export_patch_state}}) into a
##' \code{Parameters} object so that the next \code{\link{run_scm}} starts from
##' that state instead of an empty patch. The state is carried on the
##' \code{Parameters} object (rather than passed separately) so the run stays
##' self-describing and reproducible, and so the seeding survives the reset at
##' the start of every run / schedule refinement.
##'
##' @title Seed Parameters with an initial patch state
##' @param p A \code{Parameters} object (its strategies must match the exported
##'   state).
##' @param state An exported state list from \code{\link{export_patch_state}}.
##' @return The modified \code{Parameters} object.
##' @seealso \code{\link{export_patch_state}}, \code{\link{run_scm}}
##' @export
set_initial_state <- function(p, state) {
  n_spp <- length(p$strategies)
  if (length(state$n) != n_spp) {
    stop("State has ", length(state$n),
         " species but Parameters has ", n_spp, " strategies")
  }
  p$initial_state <- state$ode_state
  p$n_initial_cohorts <- as.integer(state$n)
  p$initial_node_times <- unlist(state$node_times, use.names = FALSE)
  p$initial_patch_density <- unlist(state$patch_density, use.names = FALSE)
  p$initial_pr_patch_survival <- unlist(state$pr_patch_survival, use.names = FALSE)
  p$initial_time <- state$time
  if (!is.null(state$node_schedule_times)) {
    p$node_schedule_times <- state$node_schedule_times
  }
  p
}

