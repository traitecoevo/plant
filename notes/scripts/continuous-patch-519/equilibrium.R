## The resident's equilibrium birth rate on one pulse series (#519): the birth
## rate b* at which its window-averaged fitness is one. In production regnans
## solves this; here it is a script-level root-find, to see the landscape at
## equilibrium rather than at birth rate 1.
##
## Every evaluation is one run from empty to T2 + L with the arrival window
## fixed at [T1, T2], so the realisation of the pulse series is the same at
## every birth rate and R(b) carries no window-to-window noise. Whether the
## stand had settled by T1 is reported from block means of the snapshots.
##
##   PLANT_519_SEEDS      comma-separated seeds of the series (1,2)
##   PLANT_519_T1         window start (200)
##   PLANT_519_WINDOW     window length (100)
##   PLANT_519_LIFESPAN   years past the window (50)
##   PLANT_519_BTOL       tolerance on log b (0.005)
##   PLANT_519_OUT        output directory

suppressPackageStartupMessages(library(plant))
num_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(v)
}
vec_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(strsplit(v, ",")[[1]])
}
seeds      <- vec_env("PLANT_519_SEEDS", c(1, 2))
T1         <- num_env("PLANT_519_T1", 200)
window     <- num_env("PLANT_519_WINDOW", 100)
lifespan   <- num_env("PLANT_519_LIFESPAN", 50)
btol       <- num_env("PLANT_519_BTOL", 0.005)
retire     <- num_env("PLANT_519_RETIRE", 1e-9)
lambda     <- num_env("PLANT_519_LAMBDA", 20)
mean_depth <- num_env("PLANT_519_DEPTH", 0.014)
series_tmax <- num_env("PLANT_519_SERIES_TMAX", 1000)
out_dir    <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep = " ")
  flush(stdout())
}
T2 <- T1 + window
T_end <- T2 + lifespan
resident_lma <- 0.0825
ctrl <- control(node_density_in_birth_date = TRUE, ode_step_size_max = 0.05,
                retire_survival = retire)
make_env <- function() {
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 0)
  env
}

evaluations <- list()
equilibria <- list()
for (seed in seeds) {
  set.seed(seed)
  n_pulses    <- rpois(1, lambda * series_tmax)
  pulse_time  <- sort(runif(n_pulses, 0, series_tmax))
  pulse_depth <- rexp(n_pulses, rate = 1 / mean_depth)
  k <- pulse_time > 0 & pulse_time <= T_end
  pulses <- rainfall_pulse(pulse_time[k], pulse_depth[k])

  ## R at birth rate b, and the stand's block means before the window.
  resident_R <- function(b) {
    p <- scm_base_parameters("TF24")
    p$max_patch_lifetime <- T_end
    p <- add_strategies(p, trait_matrix(resident_lma, "lma"),
                        hyperpar = TF24_hyperpar, birth_rate = b)
    p <- window_regime(p, T1, T2)
    p$node_schedule_times <- list(node_schedule_times_periodic(0, T_end, 1))
    scm <- plant:::SCM("TF24", "TF24_Env")(p, make_env(), events(node_introductions(p), pulses), ctrl)
    scm$snapshot_interval <- 1
    w <- system.time(scm$run())[["elapsed"]]
    sn <- tidy_snapshots(scm)$species
    blocks <- tapply(sn$density, pmin(floor(sn$time / 50), T1 / 50), mean)
    pre <- blocks[as.numeric(names(blocks)) < T1 / 50]
    settle <- if (length(pre) >= 2) abs(pre[length(pre)] / pre[length(pre) - 1] - 1) else NA
    R <- scm$net_reproduction_ratios
    evaluations[[length(evaluations) + 1]] <<- data.frame(
      seed = seed, b = b, R = R, wall_s = w,
      density_before_window = unname(pre[length(pre)]),
      density_in_window = mean(sn$density[sn$time >= T1 & sn$time <= T2]),
      last_block_change = unname(settle))
    stamp("seed %g: b = %.6g -> R = %.6g (%.0f s); density %.3g before the window, %.3g in it; last pre-window block moved %.1f%%",
          seed, b, R, w, pre[length(pre)], mean(sn$density[sn$time >= T1 & sn$time <= T2]), 100 * settle)
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
  f <- function(x) log(resident_R(exp(x)))
  root <- uniroot(f, c(lo, hi), f.lower = log(R_lo), f.upper = log(R_hi), tol = btol)
  b_star <- exp(root$root)
  R_star <- resident_R(b_star)
  stamp("seed %g: b* = %.6g, R(b*) = %.6g after %d root-find evaluations", seed, b_star, R_star, root$iter)
  ## How smooth R(b) is near the root: a percent either side.
  for (m in c(0.99, 1.01)) resident_R(b_star * m)
  equilibria[[length(equilibria) + 1]] <- data.frame(seed = seed, b_star = b_star, R_star = R_star)
}
e <- do.call(rbind, equilibria)
ev <- do.call(rbind, evaluations)
write.csv(e, file.path(out_dir, "equilibrium.csv"), row.names = FALSE)
write.csv(ev, file.path(out_dir, "equilibrium-evaluations.csv"), row.names = FALSE)
print(e)
stamp("done")
