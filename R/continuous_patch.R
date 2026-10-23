## Driving a continuous patch from R (plant issue #519): the thinned record a
## long run keeps, and a burn-in that stops when the stand's block means settle.

##' Flatten the thinned record of a run. With \code{scm$snapshot_interval} set
##' before the run, the SCM keeps one summary per interval instead of a copy
##' of the whole patch per leg: per species the node counts, the tallest cohort
##' carrying density, the carried density integrated over birth date and every
##' node's lifetime offspring output, plus the environment's own state.
##'
##' Read on the birth-date coordinate, which a continuous patch runs on:
##' \code{density} is the carried density weighted by each node's establishment
##' weight, a number of individuals per unit area, without the boundary node's
##' share that a census adds; on the height coordinate it is the trapezium over
##' birth dates of a density in height, which is not a count. \code{n_retired} counts the nodes carrying no density, frozen
##' or not, so it is \code{n_nodes - n_live}. The per-node \code{fecundity}
##' vectors grow with the node list, so their memory grows with the square of
##' the run's length (about 6 MB at a snapshot a year over 1000 years of one
##' species); a coarser interval cuts it.
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
  if (length(snaps) == 0L) {
    stop("this SCM kept no snapshots: set scm$snapshot_interval before its run")
  }
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
##' statistics, and the run is converged once the last \code{window} block
##' means of every statistic agree: their range is below \code{tol} times the
##' largest of them in magnitude, or below \code{tol_abs}. Judging the range
##' over a window rather than the change between neighbouring blocks is what
##' lets a steady drift accumulate into the test instead of passing block by
##' block. It stops regardless at \code{t_max}, or when the schedule is
##' consumed, and says which it was; a last block shorter than \code{block}
##' runs the patch on to \code{t_max} without being tested. Under a stochastic
##' driver the block has to span several of the driver's correlation times, or
##' the test measures the driver and not the stand. A statistic that is
##' \code{NA} in a block -- the tallest cohort of a stand with none alive --
##' cannot settle, and the function warns which.
##'
##' A settled stand statistic is not a settled fitness. A drift slower than
##' \code{tol} per \code{window * block} years passes, and an invader's fitness
##' can be far more sensitive to the stand than the statistic is: a TF24 stand
##' under constant rain at its equilibrium birth rate was still gaining 0.1 \%
##' density a century at year 5000, and the resident's fitness read off it at
##' years 300 and 500 differed sixfold. Check the quantity you will use at two
##' times.
##'
##' @title Run until the stand's block means settle
##' @param scm An \code{SCM} built and not yet run to the end, as
##'   \code{run_scm(p, ..., run = FALSE)} returns one, whose schedule reaches
##'   \code{t_max}. Its \code{snapshot_interval} is set to
##'   \code{snapshot_interval} if it is zero.
##' @param block Years per block.
##' @param tol Relative range of the last \code{window} block means below
##'   which a statistic counts as settled.
##' @param t_max Time to give up at; defaults to the run's
##'   \code{max_patch_lifetime}.
##' @param statistics Which of \code{"density"} (the species' carried density
##'   summed), \code{"height_max"} (the tallest cohort with density over the
##'   species), \code{"n_live"} (cohorts with density, summed) and
##'   \code{"environment"} (entries of the environment's state) to judge. The
##'   default leaves out \code{"height_max"} when the absorbing top layer is on,
##'   since the layer pins it at the absorbing height whatever the stand does.
##' @param snapshot_interval Years between snapshots when the SCM has none set.
##' @param environment_entries Which entries of the environment's state the
##'   \code{"environment"} statistic reads. \code{NULL} (the default) reads
##'   every entry but the cumulative flux accumulators (TF24's last five), which
##'   grow without bound, so that their relative change shrinks with time
##'   whether or not anything has settled.
##' @param window Number of consecutive block means that must agree.
##' @param tol_abs Absolute range below which a statistic counts as settled
##'   whatever its size: for one that settles at or near zero. Zero, the
##'   default, judges every statistic relatively.
##' @return A list: \code{converged}, \code{time} (where the run stands), and
##'   \code{trace}, a data frame with one row per block of the block means and
##'   the relative range of each over the window ending there (\code{rel_*},
##'   \code{NA} until \code{window} blocks have run).
##' @seealso \code{\link{tidy_snapshots}}, \code{\link{export_patch_state}}
##' @export
run_until_converged <- function(scm, block = 20, tol = 0.02,
                                t_max = scm$parameters$max_patch_lifetime,
                                statistics = c("density", "height_max", "environment"),
                                snapshot_interval = 1, environment_entries = NULL,
                                window = 3, tol_abs = 0) {
  defaulted <- missing(statistics)
  statistics <- match.arg(statistics, c("density", "height_max", "n_live", "environment"),
                          several.ok = TRUE)
  layer_on <- any(vapply(scm$patch$species, function(s) isTRUE(s$top_layer$on),
                         logical(1)))
  if (layer_on && "height_max" %in% statistics) {
    if (defaulted) {
      statistics <- setdiff(statistics, "height_max")
    } else {
      warning("with the absorbing top layer on, \"height_max\" is pinned at the ",
              "absorbing height once the layer fills, and settles whatever the ",
              "stand does", call. = FALSE)
    }
  }
  if (!(length(window) == 1L && is.finite(window) && window >= 2 &&
        window == round(window))) {
    stop("`window` must be a whole number of blocks, at least 2")
  }
  if (!(length(tol_abs) == 1L && is.finite(tol_abs) && tol_abs >= 0)) {
    stop("`tol_abs` must be a non-negative number")
  }
  if (isTRUE(scm$complete)) {
    stop("this SCM's run has already finished; build one with ",
         "run_scm(p, ..., run = FALSE) and let run_until_converged() run it")
  }
  if (!(length(block) == 1L && is.finite(block) && block > 0)) {
    stop("`block` must be a positive number of years")
  }
  if (!(length(tol) == 1L && is.finite(tol) && tol > 0)) {
    stop("`tol` must be a positive number")
  }
  if (!(length(t_max) == 1L && is.finite(t_max) && t_max > scm$time)) {
    stop("`t_max` must lie after the run's current time, ", scm$time)
  }
  n_env <- scm$patch$environment$ode_size
  if ("environment" %in% statistics && n_env == 0) {
    # In the default set it is simply not there to judge; asked for, it is a
    # request that cannot be met.
    if (defaulted) {
      statistics <- setdiff(statistics, "environment")
    } else {
      stop("this model's environment carries no state, so \"environment\" ",
           "cannot be judged; leave it out of `statistics`")
    }
  }
  if (!is.null(environment_entries)) {
    if (!"environment" %in% statistics) {
      stop("`environment_entries` is read only by the \"environment\" statistic")
    }
    if (!all(environment_entries %in% seq_len(n_env))) {
      stop("`environment_entries` must index the environment's ", n_env,
           " state entries")
    }
  } else if ("environment" %in% statistics) {
    # The accumulators are the last entries of the environment's state.
    n_acc <- length(scm$patch$state$env$soil_moist_cumulative_flux)
    environment_entries <- seq_len(n_env - n_acc)
  }
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
    if ("environment" %in% statistics) {
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
  means <- list()
  converged <- FALSE
  undefined <- character(0)
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
    undefined <- union(undefined, names(m)[is.na(m)])
    means[[k]] <- m
    rel <- rep(NA_real_, length(m))
    settled <- FALSE
    if (k >= window) {
      # One row per statistic, one column per block of the window. The range
      # is scaled by the largest magnitude in the window, not by the previous
      # block alone, so a statistic near zero is judged on its own scale and
      # one exactly zero throughout has settled.
      w <- do.call(cbind, means[(k - window + 1L):k])
      range <- apply(w, 1, function(x) max(x) - min(x))
      scale <- apply(abs(w), 1, max)
      rel <- ifelse(range == 0, 0, range / scale)
      settled <- isTRUE(all(rel < tol | range < tol_abs))
    }
    names(rel) <- paste0("rel_", names(m))
    trace[[k]] <- c(block = k, time = scm$time, m, rel)
    if (settled) {
      converged <- TRUE
      break
    }
    if (scm$complete) {
      break
    }
  }
  # What is left of the run past the last whole block.
  if (!converged && !isTRUE(scm$complete) && scm$time < t_max) {
    scm$run_to(t_max)
  }
  if (length(undefined) > 0) {
    warning("these statistics were NA in at least one block and could not ",
            "settle: ", paste(undefined, collapse = ", "),
            " (the tallest cohort of a stand with none alive is undefined)",
            call. = FALSE)
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
##' invader's fitness. The resident's own strategy scores a resident cohort's
##' lifetime output per seed, which is one at the resident's demographic
##' equilibrium (it is not \code{scm$net_reproduction_ratios}, which a
##' continuous patch's ever-growing node list does not make a per-cohort
##' quantity).
##'
##' The field is taken from the state \code{scm} has reached at the first call,
##' and held for every later one; run the resident to its fixed point first, for
##' instance with \code{\link{run_until_converged}} under constant drivers and
##' the absorbing top layer (\code{Control$absorb_epsilon}). The invaders run
##' together, as one community, since nothing they do moves a held field.
##' \code{scm} is consumed: afterwards it holds the invaders, and the field is
##' released, so a later run of it builds its own.
##'
##' The introduction grid has to make the stand a fixed point. Between two
##' introductions the field swings, and a grid that is too coarse makes the
##' stand a periodic orbit whose field, held at an introduction, is the top of
##' the swing: FF16 at its equilibrium birth rate swings by 1.8 \% in held
##' fitness at a cohort a year and 0.02 \% at four. \code{SCM$run_to()} stops
##' only at introduction instants, so to hold a field mid-interval, and check
##' the grid, end the resident's run there with \code{max_patch_lifetime}.
##'
##' @title Invader fitness at a fixed point
##' @param scm The resident's \code{SCM}, run to its fixed point.
##' @param p The resident's \code{Parameters}; each invader's strategy is built
##'   from them as \code{\link{add_mutant}} builds one.
##' @param traits A trait matrix, one row per invader.
##' @param hyperpar Hyperparameter function, as for \code{add_strategies}.
##' @param lifetime Years each invader runs for. It has to outlive the cohort:
##'   the call fails unless, at the end, none of the cohort survives, or its
##'   survival since establishment is below \code{survival_tol} and it has
##'   stopped growing. A cohort still growing may be on its way to the canopy,
##'   where most of its output lies, however few of it are left. An invader
##'   that never establishes scores zero.
##' @param survival_tol The survival since establishment below which a cohort
##'   that has stopped growing counts as dead.
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
  pm$ode_times <- numeric(0)
  pm$ode_step_sizes <- numeric(0)
  pm$ode_soil_input_slopes <- list()
  pm$ode_soil_substep_ends <- list()
  run_held <- function(years) {
    pm$max_patch_lifetime <- years
    scm$run_mutant_held(pm)
    lapply(scm$patch$species, function(sp) {
      node <- sp$nodes[[1]]
      list(state = stats::setNames(node$ode_state, node$ode_names),
           rates = stats::setNames(node$ode_rates, node$ode_names),
           output = sp$net_reproduction_ratio_by_node[[1]],
           birth_date = sp$density_in_birth_date)
    })
  }
  ## Each invader's hazard an instant after birth: what establishment cost it,
  ## so survival can be judged from establishment on, and a cohort that has
  ## barely established is not taken for a dead one.
  born <- run_held(1e-8)
  run <- run_held(lifetime)
  vapply(seq_along(run), function(i) {
    hazard_at_birth <- born[[i]]$state[["mortality"]]
    ## The fraction of seeds that establish. On the birth-date coordinate the
    ## hazard starts at zero and the newest interval grows at this fraction,
    ## which the held field keeps constant. A cohort that never established has
    ## a fraction of zero, and its output is zero however long it runs.
    pr_estab <- if (born[[i]]$birth_date) {
      born[[i]]$rates[["interval_establishment"]]
    } else {
      exp(-hazard_at_birth)
    }
    established <- pr_estab > 0
    survival <- exp(-(run[[i]]$state[["mortality"]] - hazard_at_birth))
    ## Done when none survive, or when the survivors are few and have stopped
    ## growing -- at their asymptote, or stalled in the shade, where a held
    ## field keeps them. A cohort still growing can be on its way to the canopy
    ## where most of its output lies, however few of it are left.
    height <- run[[i]]$state[["height"]]
    growing <- run[[i]]$rates[["height"]] > 1e-3 * height
    done <- !established || survival == 0 || (survival < survival_tol && !growing)
    if (!done) {
      stop(sprintf(paste("Invader %d is still alive (survival since",
                         "establishment %.3g) and %s at the end of its %g",
                         "years: raise `lifetime`"),
                   i, survival, if (growing) "growing" else "not yet below `survival_tol`",
                   lifetime))
    }
    output <- run[[i]]$output
    if (run[[i]]$birth_date) {
      output <- output * pr_estab
    }
    output * pm$strategies[[i]]$pars$S_D
  }, numeric(1))
}
