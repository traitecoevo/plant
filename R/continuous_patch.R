## Driving a continuous patch from R (plant issue #519): the thinned record a
## long run keeps, and a burn-in that stops when the stand's block means settle.

##' Flatten the thinned record of a run. With \code{scm$snapshot_interval} set
##' before the run, the SCM keeps one summary per interval instead of a copy
##' of the whole patch per leg: per species the node counts, the tallest cohort
##' carrying density, the trapezium of the carried density over the birth-date
##' grid and every node's lifetime offspring output, plus the environment's own
##' state.
##'
##' @title Tidy a run's snapshots
##' @param scm An \code{SCM} that was run with a positive
##'   \code{snapshot_interval}.
##' @return A list of three: \code{species}, a data frame with one row per
##'   snapshot and species (\code{time}, \code{species}, \code{n_nodes},
##'   \code{n_live}, \code{n_retired}, \code{height_max}, \code{density},
##'   which includes what the absorbing top layer carries, and
##'   \code{top_density});
##'   \code{environment}, a data frame with one row per snapshot of \code{time}
##'   and the environment's state entries (\code{env_1}, ...); and
##'   \code{fecundity}, a list over snapshots of lists over species of each
##'   node's lifetime offspring output.
##' @seealso \code{\link{run_until_converged}}
##' @export
tidy_snapshots <- function(scm) {
  snaps <- scm$snapshots
  species <- do.call(rbind, lapply(snaps, function(s) {
    do.call(rbind, lapply(seq_along(s$species), function(i) {
      sp <- s$species[[i]]
      data.frame(time = s$time, species = i, n_nodes = sp$n_nodes,
                 n_live = sp$n_live, n_retired = sp$n_retired,
                 height_max = sp$height_max, density = sp$density,
                 top_density = sp$top_density)
    }))
  }))
  environment <- do.call(rbind, lapply(snaps, function(s) {
    env <- s$environment
    names(env) <- if (length(env) > 0) paste0("env_", seq_along(env)) else character(0)
    c(time = s$time, env)
  }))
  list(species = species,
       environment = as.data.frame(environment),
       fecundity = lapply(snaps, function(s) lapply(s$species, `[[`, "fecundity")))
}

##' Run a patch until its stand has settled, judged on block means of its
##' snapshots: the burn-in of a continuous patch (plant issue #519), paid once
##' and then exported with \code{\link{export_patch_state}}.
##'
##' The run advances one block at a time with \code{scm$run_to()}. After each
##' block the mean over that block's snapshots is taken of the chosen
##' statistics, and the run is converged when every statistic's relative change
##' from the previous block's mean has stayed below \code{tol} for two
##' consecutive blocks. It stops regardless at \code{t_max}, or when the
##' schedule is consumed, and says which it was. Under a stochastic driver the
##' block has to span several of the driver's correlation times, or the test
##' measures the driver and not the stand.
##'
##' @title Run until the stand's block means settle
##' @param scm An \code{SCM} whose schedule reaches \code{t_max}. Its
##'   \code{snapshot_interval} is set to \code{snapshot_interval} if it is zero.
##' @param block Years per block.
##' @param tol Relative change between consecutive block means below which a
##'   statistic counts as settled.
##' @param t_max Time to give up at; defaults to the run's
##'   \code{max_patch_lifetime}.
##' @param statistics Which of \code{"density"} (the species' carried density
##'   summed), \code{"height_max"} (the tallest cohort with density over the
##'   species), \code{"n_live"} (cohorts with density, summed) and
##'   \code{"environment"} (every entry of the environment's state) to judge.
##' @param snapshot_interval Years between snapshots when the SCM has none set.
##' @param environment_entries Which entries of the environment's state the
##'   \code{"environment"} statistic reads; \code{NULL} (the default) reads them
##'   all. TF24's state is its soil water by layer followed by five cumulative
##'   flux accumulators, which grow without bound and never settle, so a TF24
##'   burn-in passes the layer indices, \code{seq_len(n_layers)}.
##' @return A list: \code{converged}, \code{time} (where the run stands), and
##'   \code{trace}, a data frame with one row per block of the block means and
##'   their relative changes (\code{rel_*}).
##' @seealso \code{\link{tidy_snapshots}}, \code{\link{export_patch_state}}
##' @export
run_until_converged <- function(scm, block = 20, tol = 0.02,
                                t_max = scm$parameters$max_patch_lifetime,
                                statistics = c("density", "height_max", "environment"),
                                snapshot_interval = 1, environment_entries = NULL) {
  statistics <- match.arg(statistics, c("density", "height_max", "n_live", "environment"),
                          several.ok = TRUE)
  if (scm$snapshot_interval <= 0) {
    scm$snapshot_interval <- snapshot_interval
  }
  if (scm$snapshot_interval > block) {
    stop("the snapshot interval must not exceed the block")
  }
  statistic <- function(s) {
    out <- c()
    if ("density" %in% statistics) {
      out <- c(out, density = sum(vapply(s$species, `[[`, numeric(1), "density")))
    }
    if ("height_max" %in% statistics) {
      h <- vapply(s$species, `[[`, numeric(1), "height_max")
      out <- c(out, height_max = if (all(is.na(h))) NA_real_ else max(h, na.rm = TRUE))
    }
    if ("n_live" %in% statistics) {
      out <- c(out, n_live = sum(vapply(s$species, `[[`, integer(1), "n_live")))
    }
    if ("environment" %in% statistics && length(s$environment) > 0) {
      env <- s$environment
      names(env) <- paste0("env_", seq_along(env))
      if (!is.null(environment_entries)) {
        env <- env[environment_entries]
      }
      out <- c(out, env)
    }
    out
  }
  t0 <- scm$time
  trace <- list()
  previous <- NULL
  streak <- 0L
  converged <- FALSE
  k <- 0L
  repeat {
    k <- k + 1L
    t_end <- t0 + k * block
    if (t_end > t_max + 1e-8) {
      break
    }
    scm$run_to(t_end)
    snaps <- scm$snapshots
    times <- vapply(snaps, `[[`, numeric(1), "time")
    in_block <- times > t_end - block & times <= t_end
    if (!any(in_block)) {
      stop("no snapshot fell inside the block ending at ", t_end)
    }
    m <- rowMeans(do.call(cbind, lapply(snaps[in_block], statistic)))
    rel <- if (is.null(previous)) rep(NA_real_, length(m)) else
      abs(m - previous) / pmax(abs(previous), .Machine$double.eps)
    names(rel) <- paste0("rel_", names(m))
    trace[[k]] <- c(block = k, time = scm$time, m, rel)
    streak <- if (!is.null(previous) && all(rel < tol)) streak + 1L else 0L
    previous <- m
    if (streak >= 2L) {
      converged <- TRUE
      break
    }
    if (scm$complete) {
      break
    }
  }
  list(converged = converged, time = scm$time,
       trace = as.data.frame(do.call(rbind, trace)))
}

##' Fitness of invaders at a resident's fixed point
##'
##' A resident at a fixed point makes one field, so an invader of it needs no
##' recording of the resident and no replay: it stands in that one field, held,
##' for its whole life. This runs each strategy in \code{traits} as one cohort
##' born at time zero into the field \code{scm} stands in now, until the cohort
##' is dead, and returns each one's lifetime offspring per seed. Every arrival
##' at a fixed point sees the same field, so that one cohort's output is the
##' invader's fitness. The resident's own strategy scores the resident's
##' reproduction ratio, which is one at its demographic equilibrium.
##'
##' The field is taken from the state \code{scm} has reached at the first call,
##' and held for every later one; run the resident to its fixed point first, for
##' instance with \code{\link{run_until_converged}} under constant drivers and
##' the absorbing top layer (\code{Control$absorb_epsilon}). The invaders run
##' together, as one community, since nothing they do moves a held field.
##' \code{scm} is consumed: afterwards it holds the invaders.
##'
##' @title Invader fitness at a fixed point
##' @param scm The resident's \code{SCM}, run to its fixed point.
##' @param p The resident's \code{Parameters}; each invader's strategy is built
##'   from them as \code{\link{add_mutant}} builds one.
##' @param traits A trait matrix, one row per invader.
##' @param hyperpar Hyperparameter function, as for \code{add_strategies}.
##' @param lifetime Years each invader runs for. It has to outlive the cohort:
##'   the call fails if a cohort's survival at the end is above
##'   \code{survival_tol}.
##' @param survival_tol The survival below which a cohort counts as dead.
##' @return Lifetime offspring per seed, one per row of \code{traits}.
##' @seealso \code{\link{run_until_converged}}
##' @export
fixed_point_fitness <- function(scm, p, traits, hyperpar = param_hyperpar(p),
                                lifetime = 3000, survival_tol = 1e-8) {
  pm <- add_mutant(p, traits, hyperpar = hyperpar,
                   birth_rate = as.list(rep(1, nrow(traits))))
  pm$patch_type <- "fixed"
  pm$max_patch_lifetime <- lifetime
  pm$node_schedule_times <- rep(list(0), nrow(traits))
  pm$initial_state <- numeric(0)
  pm$n_initial_cohorts <- integer(0)
  pm$initial_node_times <- numeric(0)
  pm$initial_patch_density <- numeric(0)
  pm$initial_pr_patch_survival <- numeric(0)
  pm$initial_time <- 0
  scm$run_mutant_held(pm)
  species <- scm$patch$species
  vapply(seq_along(species), function(i) {
    node <- species[[i]]$nodes[[1]]
    state <- stats::setNames(node$ode_state, node$ode_names)
    survival <- exp(-state[["mortality"]])
    if (survival > survival_tol) {
      stop(sprintf(paste("Invader %d is still alive at the end of its %g years",
                         "(survival %.3g): raise `lifetime`"), i, lifetime, survival))
    }
    species[[i]]$net_reproduction_ratio_by_node[[1]] * pm$strategies[[i]]$pars$S_D
  }, numeric(1))
}
