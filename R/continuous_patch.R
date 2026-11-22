## Driving a continuous patch from R (plant issue #519): the thinned record a
## long run keeps, and a burn-in that stops when the stand's block means settle.

##' Flatten the thinned record of a run. With \code{x$snapshot_interval} set
##' before the run, the deterministic runner keeps one summary per interval
##' instead of a copy
##' of the whole patch per leg: per species the node counts, the tallest cohort
##' carrying density, the trapezium of the carried density over the birth-date
##' grid and every node's lifetime offspring output, plus the environment's own
##' state.
##'
##' Read on the birth-date coordinate, which a continuous patch runs on:
##' \code{density} is the trapezium of the carried density over the nodes'
##' birth dates, a number of individuals per unit area, without the closing
##' interval to the boundary node that a census adds; on the height coordinate
##' it is the same trapezium over birth dates of a density in height, which is
##' not a count. \code{n_retired} counts the nodes carrying no density, frozen
##' or not, so it is \code{n_nodes - n_live}. The per-node \code{fecundity}
##' vectors grow with the node list, so their memory grows with the square of
##' the run's length (about 6 MB at a snapshot a year over 1000 years of one
##' species); a coarser interval cuts it.
##'
##' @title Tidy a run's snapshots
##' @param x A \code{DeterministicPatchRunner} or \code{plant_run} that was run
##'   with a positive
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
tidy_snapshots <- function(x) {
  x <- as_solver(x)
  snaps <- x$snapshots
  if (length(snaps) == 0L) {
    stop("this runner kept no snapshots: set its snapshot_interval before its run")
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
##' The run advances one block at a time with \code{x$run_to()}. After each
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
##' @param x A \code{plant_run} stopped before its end (\code{run(p, until = 0)}), or a \code{DeterministicPatchRunner} built and not yet run to the end, as
##'   \code{run(p, ..., until = 0)$solver} is, whose schedule reaches
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
##' @param snapshot_interval Years between snapshots when the deterministic
##'   runner has none set.
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
##'   \code{NA} until \code{window} blocks have run). Given a \code{plant_run},
##'   also \code{run}, the result where the run now stands.
##' @seealso \code{\link{tidy_snapshots}}, \code{\link{export_patch_state}}
##' @export
run_until_converged <- function(x, block = 20, tol = 0.02,
                                t_max = x$parameters$max_patch_lifetime,
                                statistics = c("density", "height_max", "environment"),
                                snapshot_interval = 1, environment_entries = NULL,
                                window = 3, tol_abs = 0) {
  given <- x
  x <- as_solver(x)
  checked <- check_convergence_args(x, block, tol, t_max, statistics,
                                    missing(statistics), snapshot_interval,
                                    environment_entries, window, tol_abs)
  statistic <- convergence_statistic(checked$statistics,
                                     checked$environment_entries)
  t0 <- x$time
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
    x$run_to(t_end)
    snaps <- x$snapshots
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
    trace[[k]] <- c(block = k, time = x$time, m, rel)
    if (settled) {
      converged <- TRUE
      break
    }
    if (x$complete) {
      break
    }
  }
  # What is left of the run past the last whole block.
  if (!converged && !isTRUE(x$complete) && x$time < t_max) {
    x$run_to(t_max)
  }
  if (length(undefined) > 0) {
    warning("these statistics were NA in at least one block and could not ",
            "settle: ", paste(undefined, collapse = ", "),
            " (the tallest cohort of a stand with none alive is undefined)",
            call. = FALSE)
  }
  out <- list(converged = converged, time = x$time,
              trace = as.data.frame(do.call(rbind, trace)))
  if (inherits(given, "plant_run")) {
    out$run <- new_plant_run(x, given$scenario, given$record, given$mutants)
  }
  out
}

# The arguments of run_until_converged(), checked; returns the statistics and
# environment entries it settles on. A defaulted statistic the stand cannot
# judge is dropped quietly, and an explicitly requested one is refused or warned
# about. Sets the deterministic runner's snapshot interval when it has none.
check_convergence_args <- function(solver, block, tol, t_max, statistics, defaulted,
                                   snapshot_interval, environment_entries,
                                   window, tol_abs) {
  statistics <- match.arg(statistics, c("density", "height_max", "n_live", "environment"),
                          several.ok = TRUE)
  layer_on <- any(vapply(solver$patch$species, function(s) isTRUE(s$top_layer$on),
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
    stop("`window` must be a whole number of blocks, at least 2", call. = FALSE)
  }
  if (!(length(tol_abs) == 1L && is.finite(tol_abs) && tol_abs >= 0)) {
    stop("`tol_abs` must be a non-negative number", call. = FALSE)
  }
  if (isTRUE(solver$complete)) {
    stop("this run has already finished; build one with run(p, until = 0) ",
         "and let run_until_converged() run it", call. = FALSE)
  }
  if (!(length(block) == 1L && is.finite(block) && block > 0)) {
    stop("`block` must be a positive number of years", call. = FALSE)
  }
  if (!(length(tol) == 1L && is.finite(tol) && tol > 0)) {
    stop("`tol` must be a positive number", call. = FALSE)
  }
  if (!(length(t_max) == 1L && is.finite(t_max) && t_max > solver$time)) {
    stop("`t_max` must lie after the run's current time, ", solver$time, call. = FALSE)
  }
  n_env <- solver$patch$environment$ode_size
  if ("environment" %in% statistics && n_env == 0) {
    # In the default set it is simply not there to judge; asked for, it is a
    # request that cannot be met.
    if (defaulted) {
      statistics <- setdiff(statistics, "environment")
    } else {
      stop("this model's environment carries no state, so \"environment\" ",
           "cannot be judged; leave it out of `statistics`", call. = FALSE)
    }
  }
  if (!is.null(environment_entries)) {
    if (!"environment" %in% statistics) {
      stop("`environment_entries` is read only by the \"environment\" statistic", call. = FALSE)
    }
    if (!all(environment_entries %in% seq_len(n_env))) {
      stop("`environment_entries` must index the environment's ", n_env,
           " state entries", call. = FALSE)
    }
  } else if ("environment" %in% statistics) {
    # The accumulators are the last entries of the environment's state.
    n_acc <- length(solver$patch$state$env$soil_moist_cumulative_flux)
    environment_entries <- seq_len(n_env - n_acc)
  }
  if (solver$snapshot_interval <= 0) {
    solver$snapshot_interval <- snapshot_interval
  }
  if (solver$snapshot_interval > block) {
    stop("the snapshot interval must not exceed the block", call. = FALSE)
  }
  list(statistics = statistics, environment_entries = environment_entries)
}

# The function run_until_converged() reads from each snapshot: the named
# statistics, as one numeric vector.
convergence_statistic <- function(statistics, environment_entries) {
  function(s) {
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
}

##' Fitness of invaders at a resident's fixed point
##'
##' A resident at a fixed point makes one field, so an invader of it needs no
##' recording of the resident and no pinned run: it stands in that one field,
##' constant, for its whole life. This runs each strategy in \code{traits} as one cohort
##' born at time zero into the field \code{x} stands in now, until the cohort
##' is dead, and returns each one's lifetime offspring per seed. Every arrival
##' at a fixed point sees the same field, so that one cohort's output is the
##' invader's fitness. The resident's own strategy scores a resident cohort's
##' lifetime output per seed, which is one at the resident's demographic
##' equilibrium (it is not the run's \code{net_reproduction_ratios}, which a
##' continuous patch's ever-growing node list does not make a per-cohort
##' quantity).
##'
##' The field is taken from the state \code{x} has reached at the first call,
##' and used for every later one; run the resident to its fixed point first, for
##' instance with \code{\link{run_until_converged}} under constant drivers and
##' the absorbing top layer (\code{Control$absorb_epsilon}). The invaders run
##' together, as one community, since nothing they do moves a constant field.
##' \code{x} is consumed: afterwards it holds the invaders, and the field is
##' cleared, so a later run of it builds its own. A \code{plant_run} given here
##' is refused afterwards, so pass every invader in one call.
##'
##' The introduction grid has to make the stand a fixed point. Between two
##' introductions the field swings, and a grid that is too coarse makes the
##' stand a periodic orbit whose field, taken at an introduction, is the top of
##' the swing: FF16 at its equilibrium birth rate swings by 1.8 \% in
##' constant-field fitness at a cohort a year and 0.02 \% at four. So each call
##' checks: it runs a copy of the resident on to the middle of the introduction
##' interval it stands at, scores the invaders again in the field there, and
##' warns when any invader's fitness moves by more than \code{grid_tol}
##' (relative). A warning means the stand is not at a fixed point between
##' introductions: introduce cohorts more often, or run the resident on until it
##' settles. The check costs the copy's half interval and a second run of the
##' invaders; \code{grid_tol = NULL} skips it. The copy resumes the resident, so
##' it needs the resident stopped at an introduction, as
##' \code{run(res, until = t)} and \code{run_until_converged} leave it; a
##' resident run to its end, or stopped between introductions, is warned of as
##' unchecked.
##'
##' @title Invader fitness at a fixed point
##' @param x The resident's \code{plant_run} or \code{DeterministicPatchRunner},
##'   run to its fixed point.
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
##' @param grid_tol The largest relative change in an invader's fitness between
##'   the field \code{x} stands in and the field half an introduction interval
##'   later that passes without a warning; \code{NULL} skips the check.
##' @return Lifetime offspring per seed, one per row of \code{traits}.
##' @seealso \code{\link{run_until_converged}}
##' @export
fixed_point_fitness <- function(x, p, traits, hyperpar = param_hyperpar(p),
                                lifetime = 3000, survival_tol = 1e-8,
                                grid_tol = 1e-3) {
  x <- as_solver(x)
  ## The copy is taken before the invaders' run, which consumes the resident.
  mid <- if (!is.null(grid_tol)) mid_interval_scenario(x)
  if (!is.null(grid_tol) && is.null(mid)) {
    warning("The introduction grid cannot be checked: that needs the resident ",
            "stopped at an introduction, as run(res, until = t) and ",
            "run_until_converged() leave it. Stop it there, or set ",
            "`grid_tol = NULL`", call. = FALSE)
  }
  R <- constant_field_fitness(x, p, traits, hyperpar, lifetime, survival_tol)
  if (!is.null(mid)) {
    R_mid <- constant_field_fitness(run(mid$scenario, from = mid$state)$solver,
                                    p, traits, hyperpar, lifetime, survival_tol)
    scale <- pmax(abs(R), abs(R_mid))
    swing <- ifelse(scale > 0, abs(R_mid - R) / scale, 0)
    if (any(swing > grid_tol)) {
      warning(sprintf(paste(
        "The stand is not at a fixed point between introductions: invader",
        "fitness moves by up to %.2g (relative) from the field at t = %g to the",
        "field at t = %g, mid-interval, above `grid_tol` = %g. Introduce cohorts",
        "more often, or run the resident on until it settles"),
        max(swing), mid$state$time, mid$scenario$p$max_patch_lifetime, grid_tol),
        call. = FALSE)
    }
  }
  R
}

## A copy of the resident, to run from where it stands to the middle of the
## introduction interval there: `scenario` to run, and `state` to run it from.
## run_to() stops only at introductions, so the copy ends there instead. A
## resumed run integrates only from introductions it has still to make, so the
## copy needs the resident stopped on arrival at one, as run_to() leaves it;
## NULL otherwise. The interval is the one that introduction opens, or, at the
## last, the one before it.
mid_interval_scenario <- function(solver) {
  tol <- 1e-8
  t0 <- solver$time
  p <- solver$parameters
  state <- export_patch_state(solver)
  if (!any(abs(unlist(state$node_schedule_times) - t0) <= tol)) {
    return(NULL)
  }
  tt <- sort(unique(unlist(p$node_schedule_times)))
  later <- tt[tt > t0 + tol]
  earlier <- tt[tt < t0 - tol]
  interval <- if (length(later) > 0L) {
    later[[1]] - t0
  } else if (length(earlier) > 0L) {
    t0 - earlier[[length(earlier)]]
  } else {
    return(NULL)
  }
  t_mid <- t0 + interval / 2
  ev <- state$events
  keep <- ev$time < t_mid - tol
  state$events <- Events(time = ev$time[keep], type = ev$type[keep],
                         target = ev$target[keep],
                         target_index = ev$target_index[keep],
                         params = ev$params[keep])
  p$max_patch_lifetime <- t_mid
  p$node_schedule_times <- lapply(p$node_schedule_times, function(s) s[s < t_mid])
  ## The resident's own steps are carried on its parameters; the copy takes
  ## its own.
  p$ode_times <- numeric(0)
  p$ode_step_sizes <- numeric(0)
  p$ode_soil_input_slopes <- list()
  p$ode_soil_substep_ends <- list()
  list(scenario = scenario(p, env = solver$patch$environment,
                           control = solver$control),
       state = state)
}

## Lifetime offspring per seed of each invader in `traits`, born into the field
## `x` stands in now; fixed_point_fitness() without its check.
constant_field_fitness <- function(x, p, traits, hyperpar, lifetime, survival_tol) {
  pm <- add_mutant(p, traits, hyperpar = hyperpar,
                   birth_rate = as.list(rep(1, nrow(traits))))
  pm$fitness_weighting <- "none"
  pm$max_patch_lifetime <- lifetime
  pm$node_schedule_times <- rep(list(0), nrow(traits))
  pm$initial_state <- numeric(0)
  pm$n_initial_cohorts <- integer(0)
  pm$initial_node_times <- numeric(0)
  pm$initial_patch_density <- numeric(0)
  pm$initial_pr_patch_survival <- numeric(0)
  pm$initial_time <- 0
  pm$initial_established <- FALSE
  pm$ode_times <- numeric(0)
  pm$ode_step_sizes <- numeric(0)
  pm$ode_soil_input_slopes <- list()
  pm$ode_soil_substep_ends <- list()
  run_held <- function(years) {
    pm$max_patch_lifetime <- years
    x$run_mutant_constant(pm)
    lapply(x$patch$species, function(sp) {
      node <- sp$nodes[[1]]
      state <- stats::setNames(node$ode_state, node$ode_names)
      ## On the birth-date coordinate the cohort is born per seed that
      ## establishes, and its species integrates establishment minus one over
      ## the interval after it: the first state after the nodes' blocks. The
      ## output reads it through the factor its density does, its share of that
      ## interval so far over its cell; a lone cohort's interval runs to the end
      ## of the run.
      birth_date <- sp$density_in_birth_date
      list(state = state,
           rates = stats::setNames(node$ode_rates, node$ode_names),
           output = sp$net_reproduction_ratio_by_node[[1]],
           birth_date = birth_date,
           establishment = if (birth_date)
             1 + sp$ode_rates[[sp$size * length(state) + 1L]] else NA_real_,
           factor = if (birth_date)
             exp(node$log_density - state[["log_density"]]) else 1)
    })
  }
  ## Each invader an instant after birth: what establishment cost it, so
  ## survival can be judged from establishment on, and a cohort that has
  ## barely established is not taken for a dead one.
  born <- run_held(1e-8)
  run <- run_held(lifetime)
  vapply(seq_along(run), function(i) {
    hazard_at_birth <- born[[i]]$state[["mortality"]]
    ## The fraction of seeds arriving at the birth date that establish. On the
    ## height coordinate it was paid as hazard at birth; on the birth-date one
    ## it is the rate the open interval integrates there. A cohort that never
    ## established has a fraction of zero, and its output is zero however long
    ## it runs.
    pr_estab <- if (born[[i]]$birth_date) born[[i]]$establishment else
      exp(-hazard_at_birth)
    established <- pr_estab > 0
    survival <- exp(-(run[[i]]$state[["mortality"]] - hazard_at_birth))
    ## Done when none survive, or when the survivors are few and have stopped
    ## growing -- at their asymptote, or stalled in the shade, where a constant
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
    ## Per seed arriving at the birth date: on the birth-date coordinate the
    ## output carries the mean establishment over the whole run, which a field
    ## that does not hold still (establishment decaying with patch age, say)
    ## moves away from the fraction at birth.
    output <- run[[i]]$output
    if (run[[i]]$birth_date) {
      output <- if (established && run[[i]]$factor > 0)
        output / run[[i]]$factor * pr_estab else 0
    }
    output * pm$strategies[[i]]$pars$S_D
  }, numeric(1))
}
