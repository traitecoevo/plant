## A fitness landscape under common random numbers (#519; fitness note §4.2,
## job 5): one resident burnt in and run through an arrival window on one
## pulse series, a sweep of mutant lma values invading that run in the
## resident's frozen field, and the same sweep on a second series. The relief
## of the landscape against the difference between series is what says whether
## one realisation resolves selection.
##
##   PLANT_519_SEEDS     comma-separated seeds of the series (1,2)
##   PLANT_519_LMA       comma-separated mutant lma values
##                       (0.05,0.06,0.0825,0.1,0.125,0.15,0.2)
##   PLANT_519_BIRTH     comma-separated resident birth rate per seed (1); the
##                       equilibrium values come from equilibrium.R
##   PLANT_519_T1        window start; when set, the burn-in runs to it rather
##                       than to convergence, so every series and birth rate
##                       shares one window (0, unset)
##   plus the burn-in and window settings of workflow.R.

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
lma_grid   <- vec_env("PLANT_519_LMA", c(0.05, 0.06, 0.0825, 0.1, 0.125, 0.15, 0.2))
births     <- vec_env("PLANT_519_BIRTH", rep(1, length(seeds)))
fixed_T1   <- num_env("PLANT_519_T1", 0)
if (length(births) != length(seeds)) stop("one birth rate per seed")
block      <- num_env("PLANT_519_BLOCK", 50)
tol        <- num_env("PLANT_519_TOL", 0.15)
burn_max   <- num_env("PLANT_519_BURN_MAX", 400)
window     <- num_env("PLANT_519_WINDOW", 100)
lifespan   <- num_env("PLANT_519_LIFESPAN", 50)
absorb     <- num_env("PLANT_519_ABSORB", 0.05)
lambda     <- num_env("PLANT_519_LAMBDA", 20)
mean_depth <- num_env("PLANT_519_DEPTH", 0.014)
series_tmax <- num_env("PLANT_519_SERIES_TMAX", 1000)
out_dir    <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep = " ")
  flush(stdout())
}
resident_lma <- 0.0825
stamp("plant %s; seeds %s; mutant lma %s; resident lma %g", as.character(packageVersion("plant")),
      paste(seeds, collapse = ","), paste(lma_grid, collapse = ","), resident_lma)

ctrl <- control(node_density_in_birth_date = TRUE, ode_step_size_max = 0.05,
                absorb_epsilon = absorb)
make_env <- function() {
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 0)
  env
}
parameters_for <- function(lma, lifetime, birth = 1) {
  p <- scm_base_parameters("TF24")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- lifetime
  add_strategies(p, trait_matrix(lma, "lma"), hyperpar = TF24_hyperpar, birth_rate = birth)
}

landscape <- list()
for (si in seq_along(seeds)) {
  seed <- seeds[[si]]
  birth <- births[[si]]
  set.seed(seed)
  n_pulses    <- rpois(1, lambda * series_tmax)
  pulse_time  <- sort(runif(n_pulses, 0, series_tmax))
  pulse_depth <- rexp(n_pulses, rate = 1 / mean_depth)
  intro_times <- seq(0, series_tmax, by = 1)
  pulses_in <- function(a, b) {
    k <- pulse_time > a & pulse_time <= b
    if (any(k)) rainfall_pulse(pulse_time[k], pulse_depth[k])
  }

  ## The resident's burn-in, exported where it settled.
  p <- parameters_for(resident_lma, burn_max, birth)
  p$node_schedule_times <- list(intro_times[intro_times < burn_max])
  scm <- plant:::SCM("TF24", "TF24_Env")(p, make_env(), events(node_introductions(p), pulses_in(0, burn_max)), ctrl)
  if (fixed_T1 > 0) {
    w <- system.time(scm$run_to(fixed_T1))[["elapsed"]]
    conv <- list(converged = NA)
  } else {
    w <- system.time(conv <- run_until_converged(scm, block = block, tol = tol, t_max = burn_max,
                                                 environment_entries = 1:5))[["elapsed"]]
  }
  state <- export_patch_state(scm)
  T1 <- state$time
  T2 <- T1 + window
  T_end <- T2 + lifespan
  stamp("seed %g, birth rate %.6g: burn-in %s at %g yr (%.0f s)", seed, birth,
        if (is.na(conv$converged)) "run to the fixed window start" else if (conv$converged) "converged" else "capped", T1, w)

  ## The resident through the window, from the export.
  made <- max(state$node_times[[1]])
  schedule <- intro_times[intro_times > made & intro_times < T_end]
  p2 <- parameters_for(resident_lma, T_end, birth)
  p2 <- set_initial_state(p2, state)
  p2 <- window_regime(p2, T1, T2)
  p2$node_schedule_times <- list(schedule)
  ev <- events(node_introductions(p2), pulses_in(T1, T_end))
  scm2 <- plant:::SCM("TF24", "TF24_Env")(p2, make_env(), ev, ctrl)
  w <- system.time(scm2$run())[["elapsed"]]
  R_resident <- scm2$net_reproduction_ratios
  stamp("seed %g: resident window [%g, %g] R = %.4f (%.0f s)", seed, T1, T2, R_resident, w)

  ## Mutants as a second species a million times rarer than the resident, in
  ## a full run from the same export: the field they see is the resident's to
  ## one part in a million, and the adaptive stepper is free to shrink a step
  ## their own storage pool refuses. run_mutant()'s replay by the resident's
  ## recorded step sizes cannot, and on this pulsed stand it refuses: the
  ## invader's cohorts start in a field built at their own introduction, so
  ## their storage trajectories are not the resident's cohorts', and a
  ## recorded size overshoots the pool's empty boundary (plan, design 6).
  padded <- state
  padded$n <- c(state$n, 0L)
  for (f in c("node_times", "patch_density", "pr_patch_survival", "node_schedule_times")) {
    padded[[f]] <- c(state[[f]], list(numeric(0)))
  }
  for (lma in c(resident_lma, lma_grid)) {
    pm <- scm_base_parameters("TF24")
    pm$patch_type <- "fixed"
    pm$max_patch_lifetime <- T_end
    pm <- add_strategies(pm, trait_matrix(c(resident_lma, lma), "lma"),
                         hyperpar = TF24_hyperpar, birth_rate = list(birth, birth * 1e-6))
    pm <- set_initial_state(pm, padded)
    pm <- window_regime(pm, T1, T2)
    pm$node_schedule_times <- list(schedule, schedule)
    evm <- events(node_introductions(pm), pulses_in(T1, T_end))
    scm_m <- plant:::SCM("TF24", "TF24_Env")(pm, make_env(), evm, ctrl)
    w <- system.time(scm_m$run())[["elapsed"]]
    R <- scm_m$net_reproduction_ratios
    landscape[[length(landscape) + 1]] <- data.frame(seed = seed, birth = birth, T1 = T1, lma = lma,
                                                     resident = lma == resident_lma,
                                                     R = R[[2]], R_resident_beside = R[[1]], wall_s = w)
    stamp("seed %g: mutant lma %.4f R = %.4f, resident beside it %.4f (%.0f s)%s", seed, lma, R[[2]], R[[1]], w,
          if (lma == resident_lma) sprintf("  [resident's own strategy; alone it gave %.4f]", R_resident) else "")
  }
}
d <- do.call(rbind, landscape)
out_name <- Sys.getenv("PLANT_519_NAME", "landscape")
write.csv(d, file.path(out_dir, paste0(out_name, ".csv")), row.names = FALSE)
saveRDS(d, file.path(out_dir, paste0(out_name, ".rds")))

## Relief against noise: log R by lma per series; the between-series difference
## at each lma; and whether the two series agree on where the landscape peaks.
m <- reshape(d[!duplicated(d[, c("seed", "lma")]), c("seed", "lma", "R")],
             idvar = "lma", timevar = "seed", direction = "wide")
m <- m[order(m$lma), ]
cols <- grep("^R\\.", names(m), value = TRUE)
if (length(cols) >= 2) {
  a <- log(m[[cols[1]]]); b <- log(m[[cols[2]]])
  stamp("log R relief within series: %.2f and %.2f; between-series |difference| median %.2f, max %.2f; correlation %.3f; peak lma %.4f and %.4f",
        diff(range(a)), diff(range(b)), median(abs(a - b)), max(abs(a - b)), cor(a, b),
        m$lma[which.max(a)], m$lma[which.max(b)])
  da <- diff(a); db <- diff(b)
  stamp("adjacent-lma differences agree in sign in %d of %d steps", sum(sign(da) == sign(db)), length(da))
}
print(m, digits = 5)
stamp("done")
