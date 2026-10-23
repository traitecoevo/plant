## One long patch under pulsed rainfall (#519; the semi-arid fitness note, §4):
## the pulsed TF24 stand of common.R, burnt in, exported, and resumed through a
## fitness window, with fitness the window average over arrivals. Four tasks,
## run in this order by default (PLANT_519_TASK picks some):
##
##   workflow               burn-in until the block means settle, export,
##                          resume into a window, and the arrival-output
##                          statistics: CV, autocorrelation, the window a
##                          level needs, truncation (output/workflow-*.rds)
##   equilibrium            the resident's equilibrium birth rate b* per
##                          series, by a root-find with the window fixed, so
##                          the realisation does not move with b
##                          (output/equilibrium.csv, equilibrium-evaluations.csv)
##   landscape              a mutant lma sweep on each series at birth rate 1:
##                          common random numbers (output/landscape.csv)
##   landscape_equilibrium  the same sweep at each series' b*, read from the
##                          equilibrium task, with the window fixed at
##                          PLANT_519_T1 (output/landscape-equilibrium.csv)
##
## Mutants run as a second species a million times rarer than the resident,
## in a full run from the resident's export: the field they see is the
## resident's own at every stage, at the price of the resident's full cost per
## mutant. Exact replay of the recorded resident would cost the mutant alone,
## but a replay does not yet apply event actions such as pulses (#668).
##
##   PLANT_519_SEED       seed of the workflow's series (1)
##   PLANT_519_SEEDS      comma-separated seeds of the other tasks (1,2)
##   PLANT_519_LMA        mutant lma values (0.05,0.06,0.1,0.125,0.15,0.2; the
##                        resident's own strategy is always run as well)
##   PLANT_519_BIRTH      resident birth rate per seed for `landscape` (1)
##   PLANT_519_BLOCK      years per convergence block (50)
##   PLANT_519_TOL        relative change counted as settled (0.15)
##   PLANT_519_BURN_MAX   cap on the burn-in (400)
##   PLANT_519_T1         window start for `equilibrium` and
##                        `landscape_equilibrium` (200)
##   PLANT_519_WINDOW     window length (100)
##   PLANT_519_LIFESPAN   years run past the window (50)
##   PLANT_519_BTOL       tolerance on log b in the root-find (0.005)
##   plus the stand's settings in common.R, and PLANT_519_OUT.

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "common.R"))

seed       <- num_env("PLANT_519_SEED", 1)
seeds      <- vec_env("PLANT_519_SEEDS", c(1, 2))
lma_grid   <- vec_env("PLANT_519_LMA", c(0.05, 0.06, 0.1, 0.125, 0.15, 0.2))
births     <- vec_env("PLANT_519_BIRTH", rep(1, length(seeds)))
block      <- num_env("PLANT_519_BLOCK", 50)
tol        <- num_env("PLANT_519_TOL", 0.15)
burn_max   <- num_env("PLANT_519_BURN_MAX", 400)
fixed_T1   <- num_env("PLANT_519_T1", 200)
window     <- num_env("PLANT_519_WINDOW", 100)
lifespan   <- num_env("PLANT_519_LIFESPAN", 50)
btol       <- num_env("PLANT_519_BTOL", 0.005)
intro_times <- seq(0, pulsed$series_tmax, by = 1)
tasks <- tasks_to_run(c("workflow", "equilibrium", "landscape", "landscape_equilibrium"))
stamp("plant %s at %s; tasks %s; absorb %g", as.character(packageVersion("plant")),
      find.package("plant"), paste(tasks, collapse = ", "), pulsed$absorb)

## The resident burnt in on one series, then run through the window from its
## export. With `T1` the burn-in runs to it; without, until it settles.
burn_and_window <- function(series, birth, T1 = NULL, snapshot_interval = 0) {
  p <- pulsed_parameters(pulsed$resident, burn_max, birth)
  p$node_schedule_times <- list(intro_times[intro_times < burn_max])
  scm <- run_scm(p, pulsed_env(), pulsed_control(),
                 events = events(node_introductions(p), pulses_between(series, 0, burn_max)),
                 run = FALSE)
  w_burn <- system.time(
    conv <- if (is.null(T1)) {
      run_until_converged(scm, block = block, tol = tol, t_max = burn_max,
                          statistics = c("density", "height_max", "environment"),
                          environment_entries = 1:5)
    } else {
      scm$run_to(T1)
      list(converged = NA, trace = NULL)
    })[["elapsed"]]
  state <- export_patch_state(scm)
  T1 <- state$time
  T_end <- T1 + window + lifespan
  made <- max(state$node_times[[1]])
  schedule <- intro_times[intro_times > made & intro_times < T_end]
  p2 <- set_initial_state(pulsed_parameters(pulsed$resident, T_end, birth), state)
  p2 <- fitness_window(p2, T1, T1 + window)
  p2$node_schedule_times <- list(schedule)
  scm2 <- run_scm(p2, pulsed_env(), pulsed_control(),
                  events = events(node_introductions(p2), pulses_between(series, T1, T_end)),
                  run = FALSE)
  scm2$snapshot_interval <- snapshot_interval
  w_window <- system.time(scm2$run())[["elapsed"]]
  list(burn = scm, conv = conv, w_burn = w_burn, state = state, window = scm2,
       w_window = w_window, T1 = T1, T2 = T1 + window, T_end = T_end,
       schedule = schedule, p2 = p2)
}

## ---- workflow -------------------------------------------------------------
if ("workflow" %in% tasks) {
  series <- pulse_series_for(seed)
  r <- burn_and_window(series, birth = 1, snapshot_interval = 1)
  sp0 <- r$burn$patch$species[[1]]
  stamp("burn-in: %s at %g yr in %.1f s; %d cohorts, %d with density, %d retired",
        if (r$conv$converged) "converged" else "NOT converged", r$T1, r$w_burn, sp0$size,
        sum(!sp0$retired), sum(sp0$retired))
  print(round(r$conv$trace, 4))
  saveRDS(list(trace = r$conv$trace, converged = r$conv$converged, time = r$T1,
               wall = r$w_burn, snapshots = tidy_snapshots(r$burn), state = r$state),
          file.path(out_dir, "workflow-burnin.rds"))

  sp <- r$window$patch$species[[1]]
  S_D <- r$p2$strategies[[1]]$pars$S_D
  t <- sp$node_times
  out <- sp$net_reproduction_ratio_by_node * S_D
  inside <- t >= r$T1 & t <= r$T2
  k <- which(inside)
  stamp("window [%g, %g] run to %g in %.1f s: R = %.4f at birth rate 1 (replica %.4f); %d cohorts in the window, %d seeded from the burn-in at weight %g",
        r$T1, r$T2, r$T_end, r$w_window, r$window$net_reproduction_ratios,
        trapezium(t[k], out[k]) / (t[max(k)] - t[min(k)]), sum(inside), sum(t < r$T1),
        max(sp$patch_densities[t < r$T1]))

  ## How variable is one arrival's output along the window?
  x <- out[inside]
  t0 <- t[inside]
  n <- length(x)
  rho <- if (n > 2) cor(x[-1], x[-n]) else NA
  n_eff <- n * (1 - rho) / (1 + rho)
  cv <- sd(x) / mean(x)
  half <- function(a, b) { j <- which(t >= a & t <= b); trapezium(t[j], out[j]) / (t[max(j)] - t[min(j)]) }
  stamp("arrival output: mean %.4f, CV %.2f, zero in %d of %d years, lag-1 autocorrelation %.2f, n_eff %.0f, relative SE of the window mean %.3f",
        mean(x), cv, sum(x == 0), n, rho, n_eff, cv / sqrt(n_eff))
  stamp("half windows: [%g, %g] %.4f, [%g, %g] %.4f", r$T1, r$T1 + window / 2,
        half(r$T1, r$T1 + window / 2), r$T1 + window / 2, r$T2, half(r$T1 + window / 2, r$T2))
  stamp("window of years needed for a 1%% relative SE at this CV and autocorrelation: %.0f; for 5%%: %.0f",
        cv^2 * (1 + rho) / (1 - rho) / 0.01^2, cv^2 * (1 + rho) / (1 - rho) / 0.05^2)

  ## Truncation: how much of its final output a cohort born in the window's
  ## first decade had reached by each age, off the snapshots' per-node output.
  sn <- tidy_snapshots(r$window)
  fec <- lapply(sn$fecundity, function(f) f[[1]])
  snap_times <- unique(sn$species$time)
  n_seed <- sum(t < r$T1)
  accrued <- sapply(1:10, function(j) {
    idx <- n_seed + j
    born <- r$T1 + j - 1
    final <- fec[[length(fec)]][idx]
    sapply(c(10, 25, 50), function(age) {
      i <- which(snap_times >= born + age)[1]
      if (is.na(i) || length(fec[[i]]) < idx || final == 0) NA else fec[[i]][idx] / final
    })
  })
  stamp("median fraction of a cohort's final output reached by age 10, 25, 50: %s",
        paste(round(apply(accrued, 1, median, na.rm = TRUE), 3), collapse = ", "))
  stamp("output by decade of arrival: %s",
        paste(round(tapply(x, floor((t0 - r$T1) / 10), mean), 1), collapse = " "))
  saveRDS(list(T1 = r$T1, T2 = r$T2, T_end = r$T_end, R_window = r$window$net_reproduction_ratios,
               wall = r$w_window, t0 = t0, r = x, rho = rho, n_eff = n_eff,
               rel_se = cv / sqrt(n_eff), accrued = accrued, snapshots = sn),
          file.path(out_dir, "workflow-window.rds"))
}

## ---- equilibrium -----------------------------------------------------------
## Every evaluation is one run from empty to T2 + L with the window fixed at
## [T1, T2], so the realisation of the series is the same at every birth rate
## and R(b) carries no window-to-window noise. In production regnans solves
## this fixed point; here it is a script-level root-find.
if ("equilibrium" %in% tasks) {
  T2 <- fixed_T1 + window
  T_end <- T2 + lifespan
  evaluations <- list()
  equilibria <- list()
  for (s in seeds) {
    series <- pulse_series_for(s)
    pulses <- pulses_between(series, 0, T_end)
    resident_R <- function(b) {
      p <- pulsed_parameters(pulsed$resident, T_end, b)
      p <- fitness_window(p, fixed_T1, T2)
      p$node_schedule_times <- list(node_schedule_times_periodic(0, T_end, 1))
      scm <- run_scm(p, pulsed_env(), pulsed_control(),
                     events = events(node_introductions(p), pulses), run = FALSE)
      scm$snapshot_interval <- 1
      w <- system.time(scm$run())[["elapsed"]]
      sn <- tidy_snapshots(scm)$species
      blocks <- tapply(sn$density, pmin(floor(sn$time / 50), fixed_T1 / 50), mean)
      pre <- blocks[as.numeric(names(blocks)) < fixed_T1 / 50]
      settle <- if (length(pre) >= 2) abs(pre[length(pre)] / pre[length(pre) - 1] - 1) else NA
      R <- scm$net_reproduction_ratios
      in_window <- mean(sn$density[sn$time >= fixed_T1 & sn$time <= T2])
      evaluations[[length(evaluations) + 1]] <<- data.frame(
        seed = s, b = b, R = R, wall_s = w, density_before_window = unname(pre[length(pre)]),
        density_in_window = in_window, last_block_change = unname(settle))
      stamp("seed %g: b = %.6g -> R = %.6g (%.0f s); density %.3g before the window, %.3g in it; last pre-window block moved %.1f%%",
            s, b, R, w, pre[length(pre)], in_window, 100 * settle)
      R
    }
    ## Bracket the root of log R in log b, stepping up by a decade from b = 1.
    lo <- 0
    R_lo <- resident_R(exp(lo))
    if (!(R_lo > 1)) stop("R <= 1 at b = 1; bracket downward instead")
    hi <- lo
    R_hi <- R_lo
    while (R_hi > 1) {
      lo <- hi
      R_lo <- R_hi
      hi <- hi + log(10)
      R_hi <- resident_R(exp(hi))
      if (hi > log(1e6)) stop("no root below b = 1e6")
    }
    root <- uniroot(function(x) log(resident_R(exp(x))), c(lo, hi),
                    f.lower = log(R_lo), f.upper = log(R_hi), tol = btol)
    b_star <- exp(root$root)
    R_star <- resident_R(b_star)
    stamp("seed %g: b* = %.6g, R(b*) = %.6g after %d root-find evaluations", s, b_star, R_star, root$iter)
    ## How smooth R(b) is near the root: a percent either side.
    for (m in c(0.99, 1.01)) resident_R(b_star * m)
    equilibria[[length(equilibria) + 1]] <- data.frame(seed = s, b_star = b_star, R_star = R_star)
  }
  e <- do.call(rbind, equilibria)
  write.csv(e, file.path(out_dir, "equilibrium.csv"), row.names = FALSE)
  write.csv(do.call(rbind, evaluations), file.path(out_dir, "equilibrium-evaluations.csv"),
            row.names = FALSE)
  print(e)
}

## ---- landscape and landscape_equilibrium -------------------------------------
## Mutants as a second species beside the resident, from the resident's export.
landscape_on <- function(seeds, births, T1, name) {
  rows <- list()
  for (si in seq_along(seeds)) {
    s <- seeds[[si]]
    birth <- births[[si]]
    series <- pulse_series_for(s)
    r <- burn_and_window(series, birth, T1)
    R_resident <- r$window$net_reproduction_ratios
    stamp("seed %g, birth rate %.6g: burn-in %s at %g yr (%.0f s); resident window [%g, %g] R = %.4f (%.0f s)",
          s, birth, if (is.na(r$conv$converged)) "run to the fixed window start" else
            if (r$conv$converged) "converged" else "capped", r$T1, r$w_burn, r$T1, r$T2, R_resident,
          r$w_window)
    ## The second species' block goes between the first's and the environment's:
    ## no nodes, but its top-layer block when the layer is on (the strategy's
    ## states plus five; a node's is the states plus two).
    padded <- r$state
    padded$n <- c(r$state$n, 0L)
    for (f in c("node_times", "patch_density", "pr_patch_survival", "node_schedule_times")) {
      padded[[f]] <- c(r$state[[f]], list(numeric(0)))
    }
    n_env <- pulsed_env()$ode_size
    top_size <- if (pulsed$absorb > 0) {
      length(r$window$patch$species[[1]]$new_node$ode_names) - 2L + 5L
    } else 0L
    padded$ode_state <- c(head(r$state$ode_state, -n_env), rep(0, top_size),
                          tail(r$state$ode_state, n_env))
    for (lma in unique(c(pulsed$resident, lma_grid))) {
      pm <- scm_base_parameters("TF24")
      pm$patch_type <- "fixed"
      pm$max_patch_lifetime <- r$T_end
      pm <- add_strategies(pm, trait_matrix(c(pulsed$resident, lma), "lma"),
                           hyperpar = TF24_hyperpar, birth_rate = list(birth, birth * 1e-6))
      pm <- set_initial_state(pm, padded)
      pm <- fitness_window(pm, r$T1, r$T2)
      pm$node_schedule_times <- list(r$schedule, r$schedule)
      scm_m <- run_scm(pm, pulsed_env(), pulsed_control(),
                       events = events(node_introductions(pm), pulses_between(series, r$T1, r$T_end)),
                       run = FALSE)
      w <- system.time(scm_m$run())[["elapsed"]]
      R <- scm_m$net_reproduction_ratios
      rows[[length(rows) + 1]] <- data.frame(seed = s, birth = birth, T1 = r$T1, lma = lma,
                                             resident = lma == pulsed$resident, R = R[[2]],
                                             R_resident_beside = R[[1]],
                                             R_resident_alone = R_resident, wall_s = w)
      stamp("seed %g: mutant lma %.4f R = %.4f, resident beside it %.4f (%.0f s)%s", s, lma, R[[2]],
            R[[1]], w, if (lma == pulsed$resident) sprintf("  [resident's own strategy; alone it gave %.4f]", R_resident) else "")
    }
  }
  d <- do.call(rbind, rows)
  write.csv(d, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)

  ## Relief against noise: the ratio to the resident per series, the gap
  ## between series in it, and whether the series agree on the landscape.
  d$rel <- d$R / d$R_resident_alone
  w <- reshape(d[, c("seed", "lma", "R", "rel")], idvar = "lma", timevar = "seed",
               direction = "wide")
  w <- w[order(w$lma), ]
  rel_cols <- grep("^rel\\.", names(w), value = TRUE)
  if (length(rel_cols) >= 2) {
    w$gap_pct <- 100 * abs(w[[rel_cols[2]]] / w[[rel_cols[1]]] - 1)
    a <- log(w[[rel_cols[1]]]); b <- log(w[[rel_cols[2]]])
    g <- function(col) {
      (log(w[[col]][w$lma == 0.1]) - log(w[[col]][w$lma == 0.06])) / 0.04
    }
    stamp("d log R / d lma across [0.06, 0.10]: %.1f and %.1f; adjacent-lma differences agree in sign in %d of %d steps",
          g(rel_cols[1]), g(rel_cols[2]), sum(sign(diff(a)) == sign(diff(b))), length(a) - 1)
  }
  print(format(w, digits = 5), row.names = FALSE)
}
if ("landscape" %in% tasks) {
  landscape_on(seeds, births, T1 = NULL, name = "landscape")
}
if ("landscape_equilibrium" %in% tasks) {
  e <- read.csv(file.path(out_dir, "equilibrium.csv"))
  landscape_on(e$seed, e$b_star, T1 = fixed_T1, name = "landscape-equilibrium")
}
stamp("done")
