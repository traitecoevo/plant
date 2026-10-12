## The continuous-patch workflow end to end on the pulsed TF24 stand (#519):
## burn in until the stand's block means settle, export, resume into an
## arrival window, read the resident's window-averaged fitness, and measure
## how variable one cohort's lifetime output is along the window -- which is
## what sets the window length the fitness note's §4.5 asks for.
##
## Same stand and the same pulse series as measure.R (seed 1 over 1000 yr).
##   PLANT_519_BLOCK      years per convergence block (50)
##   PLANT_519_TOL        relative change between block means that counts as
##                        settled (0.15; the driver's fluctuation sets the floor)
##   PLANT_519_BURN_MAX   cap on the burn-in (400)
##   PLANT_519_WINDOW     window length (100)
##   PLANT_519_LIFESPAN   years run past the window so its last arrivals live
##                        out their output (50; the pilot's cohorts die by 30)
##   PLANT_519_RETIRE     control$retire_density (1e-9)
##   PLANT_519_OUT        output directory

suppressPackageStartupMessages(library(plant))
num_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(v)
}
block      <- num_env("PLANT_519_BLOCK", 50)
tol        <- num_env("PLANT_519_TOL", 0.15)
burn_max   <- num_env("PLANT_519_BURN_MAX", 400)
window     <- num_env("PLANT_519_WINDOW", 100)
lifespan   <- num_env("PLANT_519_LIFESPAN", 50)
retire     <- num_env("PLANT_519_RETIRE", 1e-9)
lambda     <- num_env("PLANT_519_LAMBDA", 20)
mean_depth <- num_env("PLANT_519_DEPTH", 0.014)
seed       <- num_env("PLANT_519_SEED", 1)
series_tmax <- num_env("PLANT_519_SERIES_TMAX", 1000)
out_dir    <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep = " ")
  flush(stdout())
}
stamp("plant %s at %s; block %g, tol %g, burn-in cap %g, window %g, lifespan %g, retire %g",
      as.character(packageVersion("plant")), find.package("plant"), block, tol, burn_max,
      window, lifespan, retire)

set.seed(seed)
n_pulses    <- rpois(1, lambda * series_tmax)
pulse_time  <- sort(runif(n_pulses, 0, series_tmax))
pulse_depth <- rexp(n_pulses, rate = 1 / mean_depth)
intro_times <- seq(0, series_tmax, by = 1)

base_parameters <- function(lifetime) {
  p <- scm_base_parameters("TF24")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- lifetime
  add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = TF24_hyperpar, birth_rate = 1)
}
ctrl <- control(node_density_in_birth_date = TRUE, ode_step_size_max = 0.05,
                retire_density = retire)
make_env <- function() {
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 0)
  env
}
pulses_in <- function(a, b) {
  k <- pulse_time > a & pulse_time <= b
  if (any(k)) rainfall_pulse(pulse_time[k], pulse_depth[k])
}

## ---- 1. burn-in, until the block means settle ------------------------------
p <- base_parameters(burn_max)
p$node_schedule_times <- list(intro_times[intro_times < burn_max])
scm <- plant:::SCM("TF24", "TF24_Env")(p, make_env(), events(node_introductions(p), pulses_in(0, burn_max)), ctrl)
w_burn <- system.time(conv <- run_until_converged(
  scm, block = block, tol = tol, t_max = burn_max,
  statistics = c("density", "height_max", "environment"), environment_entries = 1:5))[["elapsed"]]
stamp("burn-in: %s at %g yr in %.1f s; %d cohorts, %d with density, %d retired",
      if (conv$converged) "converged" else "NOT converged", conv$time, w_burn,
      scm$patch$species[[1]]$size, sum(exp(scm$patch$species[[1]]$log_densities_state) > 0),
      sum(scm$patch$species[[1]]$retired))
print(round(conv$trace, 4))
state <- export_patch_state(scm)
burn_snapshots <- tidy_snapshots(scm)
saveRDS(list(trace = conv$trace, converged = conv$converged, time = conv$time, wall = w_burn,
             snapshots = burn_snapshots, state = state),
        file.path(out_dir, "workflow-burnin.rds"))

## ---- 2. the window run, resumed from the export ----------------------------
T1 <- state$time
T2 <- T1 + window
T_end <- T2 + lifespan
p2 <- base_parameters(T_end)
p2 <- set_initial_state(p2, state)
p2 <- window_regime(p2, T1, T2)
made <- max(state$node_times[[1]])
p2$node_schedule_times <- list(intro_times[intro_times > made & intro_times < T_end])
scm2 <- plant:::SCM("TF24", "TF24_Env")(p2, make_env(), events(node_introductions(p2), pulses_in(T1, T_end)), ctrl)
scm2$snapshot_interval <- 1
w_win <- system.time(scm2$run())[["elapsed"]]
R_window <- scm2$net_reproduction_ratios
sp <- scm2$patch$species[[1]]
S_D <- p2$strategies[[1]]$pars$S_D
t <- sp$node_times
out <- sp$net_reproduction_ratio_by_node * S_D
inside <- t >= T1 & t <= T2
## The regime's number is the trapezium-weighted mean of the window's arrivals,
## and nothing of the seeded burn-in cohorts: replicated here as a check.
trap <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
d <- ifelse(inside, 1 / window, 0)
stamp("window [%g, %g] run to %g in %.1f s: R = %.4f at birth rate 1 (replica %.4f); %d cohorts in the window, %d seeded from the burn-in at weight %g",
      T1, T2, T_end, w_win, R_window, trap(t, out * d) / trap(t, d), sum(inside),
      sum(t < T1), max(sp$patch_densities[t < T1]))

## ---- 3. how variable is one arrival's output along the window? --------------
r <- out[inside]
t0 <- t[inside]
n <- length(r)
rho <- if (n > 2) cor(r[-1], r[-n]) else NA
n_eff <- n * (1 - rho) / (1 + rho)
rel_se <- sd(r) / mean(r) / sqrt(n_eff)
half <- function(a, b) { d <- ifelse(t >= a & t <= b, 1, 0); trap(t, out * d) / trap(t, d) }
stamp("arrival output: mean %.4f, CV %.2f, zero in %d of %d years, lag-1 autocorrelation %.2f, n_eff %.0f, relative SE of the window mean %.3f",
      mean(r), sd(r) / mean(r), sum(r == 0), n, rho, n_eff, rel_se)
stamp("half windows: [%g, %g] %.4f, [%g, %g] %.4f", T1, T1 + window / 2, half(T1, T1 + window / 2),
      T1 + window / 2, T2, half(T1 + window / 2, T2))
stamp("window of years needed for a 1%% relative SE at this CV and autocorrelation: %.0f; for 5%%: %.0f",
      (sd(r) / mean(r))^2 * (1 + rho) / (1 - rho) / 0.01^2, (sd(r) / mean(r))^2 * (1 + rho) / (1 - rho) / 0.05^2)
## Truncation: how much of its final output a cohort born in the window's first
## decade had reached by each age, off the per-node output in the snapshots.
sn <- tidy_snapshots(scm2)
fec <- lapply(sn$fecundity, function(f) f[[1]])
snap_times <- unique(sn$species$time)
n_seed <- sum(t < T1)
accrued <- sapply(1:10, function(j) {
  idx <- n_seed + j
  born <- T1 + j - 1
  final <- fec[[length(fec)]][idx]
  sapply(c(10, 25, 50), function(age) {
    k <- which(snap_times >= born + age)[1]
    if (is.na(k) || length(fec[[k]]) < idx || final == 0) NA else fec[[k]][idx] / final
  })
})
stamp("median fraction of a cohort's final output reached by age 10, 25, 50: %s",
      paste(round(apply(accrued, 1, median, na.rm = TRUE), 3), collapse = ", "))
stamp("output by decade of arrival: %s", paste(round(tapply(r, floor((t0 - T1) / 10), mean), 1), collapse = " "))
saveRDS(list(T1 = T1, T2 = T2, T_end = T_end, R_window = R_window, wall = w_win,
             t0 = t0, r = r, rho = rho, n_eff = n_eff, rel_se = rel_se,
             accrued = accrued, snapshots = sn),
        file.path(out_dir, "workflow-window.rds"))
stamp("done")
