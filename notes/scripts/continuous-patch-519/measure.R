## Pilot cost measurement for the continuous patch (plant #519).
##
## One TF24 patch, no disturbance regime, driven by a synthetic stationary
## rainfall-pulse series through the event queue (#632), run as long as a wall
## budget allows. Reports wall time per simulated year, accepted steps per year,
## node counts (introduced and live), the stand's height and density, soil water
## by layer and the water budget per chunk, and whether the stand approaches a
## stochastic steady state.
##
## The run is chunked through the resume path (export_patch_state() /
## set_initial_state()), which is the mechanism the design rests on, so the
## pilot also measures how far a resumed run departs from a straight one.
##
## Settings come from the environment so the same script serves a smoke test
## and the full pilot:
##   PLANT_519_BUDGET_S  wall budget in seconds (7200)
##   PLANT_519_CHUNK     chunk length in years (10)
##   PLANT_519_TMAX      horizon in years (1000)
##   PLANT_519_LAMBDA    pulses per year (20)
##   PLANT_519_DEPTH     mean pulse depth in m, exponential (0.014 -> 280 mm/yr)
##   PLANT_519_SEED      seed of the pulse series (1)
##   PLANT_519_PROBE     1 to run the resume-vs-straight probe first (1)
##   PLANT_519_ABSORB    control$absorb_epsilon, the absorbing top layer; 0 (the
##                       default) absorbs nothing
##   PLANT_519_SERIES_TMAX  horizon the pulse series is drawn over (TMAX); set it
##                       to a longer run's horizon to reuse that run's series
##   PLANT_519_OUT       output directory (notes/scripts/continuous-patch-519/output)
##
## Run from the plant root with the library that holds the pinned odelia and
## phylloptim:  R_LIBS=<lib> Rscript --no-init-file notes/scripts/continuous-patch-519/measure.R

suppressPackageStartupMessages(library(plant))

num_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(v)
}
budget_s   <- num_env("PLANT_519_BUDGET_S", 7200)
chunk      <- num_env("PLANT_519_CHUNK", 10)
t_max      <- num_env("PLANT_519_TMAX", 1000)
lambda     <- num_env("PLANT_519_LAMBDA", 20)
mean_depth <- num_env("PLANT_519_DEPTH", 0.014)
seed       <- num_env("PLANT_519_SEED", 1)
do_probe   <- num_env("PLANT_519_PROBE", 1) > 0
absorb     <- num_env("PLANT_519_ABSORB", 0)
series_tmax <- num_env("PLANT_519_SERIES_TMAX", t_max)
out_dir    <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep = " ")
  flush(stdout())
}

stamp("plant %s at %s", as.character(packageVersion("plant")), find.package("plant"))
stamp("odelia %s, phylloptim %s", as.character(packageVersion("odelia")),
      as.character(packageVersion("phylloptim")))
stamp("git %s", tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                        error = function(e) "?"))
stamp("budget %.0f s, chunk %g yr, horizon %g yr, %g pulses/yr of mean %.3f m, seed %g, absorb_epsilon %g",
      budget_s, chunk, t_max, lambda, mean_depth, seed, absorb)

## ---- the stand ---------------------------------------------------------------

## patch_type other than "meta-population" selects No_Disturbance: patch density
## 1 and survival 1, so lifetime offspring are unweighted.
base_parameters <- function() {
  p <- scm_base_parameters("TF24")
  p$patch_type <- "fixed"
  add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = TF24_hyperpar,
                 birth_rate = 1)
}

## Birth-date coordinate (#590) and the measured step cap for TF24's stiff soil
## column (#599): the cap bounds the step a leg inherits across the boundary a
## pulse makes, and pulses make many.
ctrl <- control(node_density_in_birth_date = TRUE, ode_step_size_max = 0.05,
                absorb_epsilon = absorb)

## Rainfall as pulses alone: the continuous driver is turned off, or the run
## would receive both.
make_env <- function() {
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 0)
  env
}

## ---- the forcing, baked in up front ----------------------------------------

## A stationary marked Poisson process: exponential inter-arrival times and
## exponential depths. One realisation, drawn once, shared by every segment.
set.seed(seed)
## pulse_series_poisson() draws count, times, then amounts, so a seed rebuilds
## the series this script's committed numbers were measured on.
series      <- pulse_series_poisson(0, series_tmax, rate = lambda, mean_amount = mean_depth)
pulse_time  <- series$time
pulse_depth <- series$amount
n_pulses    <- nrow(series)
stamp("%d pulses over %g yr, %.1f mm/yr", n_pulses, series_tmax,
      1000 * sum(pulse_depth) / series_tmax)

## Periodic introduction: one cohort per year, on the integers.
intro_times <- seq(0, t_max, by = 1)

## ---- one segment --------------------------------------------------------------

## Run the patch from t0 to t1, from empty or from an exported state. The
## segment's own introductions and pulses are the ones inside (t0, t1]; an
## introduction AT t0 belongs to the segment that starts there (only the first
## segment has one, at 0; chunk ends sit on half-years so no later one does).
run_segment <- function(t0, t1, state = NULL) {
  p <- base_parameters()
  p$max_patch_lifetime <- t1
  if (!is.null(state)) {
    p <- set_initial_state(p, state)
  }
  ## After add_strategies(), which is what makes Parameters re-validate; the
  ## residual schedule set_initial_state() installed is replaced by this one.
  intro <- intro_times[intro_times >= t0 & intro_times < t1]
  p$node_schedule_times <- list(intro)
  in_seg <- pulse_time > t0 & pulse_time <= t1
  ## events() drops NULL parts; a segment with neither (a short tail) runs the
  ## empty schedule, which is a no-op, and its summary says so with zero steps.
  ev <- events(if (length(intro) > 0) node_introductions(p),
               if (any(in_seg)) rainfall_pulse(pulse_time[in_seg], pulse_depth[in_seg]))
  run_scm(p, make_env(), ctrl, events = ev)
}

## Per-node state of species 1 off the flat ODE vector: the nodes come first,
## one block of node_ode_size entries each, in the order ode_names gives.
node_matrix <- function(scm) {
  patch <- scm$patch
  sp <- patch$species[[1]]
  names <- sp$new_node$ode_names
  n <- sp$size
  k <- length(names)
  m <- matrix(patch$ode_state[seq_len(n * k)], nrow = k,
              dimnames = list(names, NULL))
  m
}

trapezium <- function(x, y) {
  if (length(x) < 2) return(0)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

summarise_segment <- function(scm, t0, t1, wall) {
  patch <- scm$patch
  sp <- patch$species[[1]]
  m <- node_matrix(scm)
  times <- sp$node_times
  h <- m["height", ]
  nu <- exp(m["log_density", ])          # density in birth date
  fec <- m["offspring_produced_survival_weighted", ]
  live <- nu > 1e-6 * max(nu)
  retired <- sp$retired
  env <- patch$environment
  theta <- env$get_soil_water_state()
  flux <- env$get_soil_water_state_cumulative_flux()
  log <- scm$event_log
  is_pulse <- log$type == "resource_pulse"
  accepted <- sum(vapply(log$applied[is_pulse], `[[`, numeric(1), 1))
  shed <- sum(vapply(log$applied[is_pulse], `[[`, numeric(1), 2))
  data.frame(
    t0 = t0, t1 = t1, years = t1 - t0, wall_s = wall,
    s_per_year = wall / (t1 - t0),
    steps = length(scm$ode_times), steps_per_year = length(scm$ode_times) / (t1 - t0),
    n_nodes = sp$size, n_live = sum(live), n_retired = sum(retired),
    height_max = max(h), height_mean_live = sum(h[live] * nu[live]) / sum(nu[live]),
    density_total = trapezium(times, nu),
    fecundity_total = trapezium(times, fec),
    offspring_production = scm$offspring_production[[1]],
    theta_1 = theta[[1]], theta_2 = theta[[2]], theta_3 = theta[[3]],
    theta_4 = theta[[4]], theta_5 = theta[[5]],
    sum_rainfall = flux[[1]], sum_infiltration = flux[[2]],
    sum_drainage = flux[[3]], sum_uptake = flux[[4]], sum_pulse_runoff = flux[[5]],
    pulses = sum(is_pulse), pulse_accepted_m = accepted, pulse_shed_m = shed
  )
}

## ---- probe: resume versus straight -------------------------------------------

if (do_probe) {
  stamp("probe: straight 0 -> 2.5 versus 0 -> 1.5 -> 2.5")
  w1 <- system.time(straight <- run_segment(0, 2.5))[["elapsed"]]
  w2 <- system.time(first <- run_segment(0, 1.5))[["elapsed"]]
  w3 <- system.time(resumed <- run_segment(1.5, 2.5, export_patch_state(first)))[["elapsed"]]
  a <- straight$patch$ode_state
  b <- resumed$patch$ode_state
  rel <- abs(a - b) / pmax(abs(a), 1e-300)
  stamp("probe: %d vs %d state entries, max |rel diff| %.3e (median %.3e), offspring %.6g vs %.6g, wall %.1f vs %.1f + %.1f s",
        length(a), length(b), max(rel), median(rel),
        straight$offspring_production[[1]], resumed$offspring_production[[1]], w1, w2, w3)
  saveRDS(list(straight = a, resumed = b, times = straight$patch$species[[1]]$node_times),
          file.path(out_dir, "probe.rds"))
}

## ---- the long run ---------------------------------------------------------------

rows <- list()
state <- NULL
t0 <- 0
k <- 0
elapsed <- 0
failure <- NULL
repeat {
  if (t0 >= t_max) break
  k <- k + 1
  t1 <- min(k * chunk + 0.5, t_max)
  res <- tryCatch({
    wall <- system.time(scm <- run_segment(t0, t1, state))[["elapsed"]]
    list(ok = TRUE, scm = scm, wall = wall)
  }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  if (!res$ok) {
    failure <- list(t0 = t0, t1 = t1, message = res$message)
    stamp("FAILED in segment %g -> %g: %s", t0, t1, res$message)
    break
  }
  row <- summarise_segment(res$scm, t0, t1, res$wall)
  rows[[k]] <- row
  table <- do.call(rbind, rows)
  write.csv(table, file.path(out_dir, "chunks.csv"), row.names = FALSE)
  stamp("%6.1f yr: %6.1f s (%5.2f s/yr, %5.0f steps/yr) nodes %4d live %4d retired %4d h_max %5.2f dens %.3g R0 %.3g theta %.3f/%.3f/%.3f/%.3f/%.3f shed %.0f%%",
        t1, res$wall, row$s_per_year, row$steps_per_year, row$n_nodes, row$n_live, row$n_retired,
        row$height_max, row$density_total, row$offspring_production,
        row$theta_1, row$theta_2, row$theta_3, row$theta_4, row$theta_5,
        100 * row$pulse_shed_m / max(row$pulse_accepted_m + row$pulse_shed_m, 1e-12))
  state <- export_patch_state(res$scm)
  saveRDS(state, file.path(out_dir, "last_state.rds"))
  saveRDS(node_matrix(res$scm), file.path(out_dir, "last_nodes.rds"))
  elapsed <- elapsed + res$wall
  ## The next chunk costs at least what this one did, and more as nodes grow.
  if (elapsed + 1.5 * res$wall > budget_s) {
    stamp("budget: %.0f s used, next chunk predicted >= %.0f s; stopping at %g yr",
          elapsed, 1.5 * res$wall, t1)
    break
  }
  t0 <- t1
}

if (length(rows) > 0) {
  table <- do.call(rbind, rows)
  saveRDS(list(chunks = table, failure = failure,
               settings = list(budget_s = budget_s, chunk = chunk, t_max = t_max,
                               lambda = lambda, mean_depth = mean_depth, seed = seed,
                               absorb_epsilon = absorb),
               versions = list(plant = as.character(packageVersion("plant")),
                               odelia = as.character(packageVersion("odelia")),
                               phylloptim = as.character(packageVersion("phylloptim")))),
          file.path(out_dir, "measurement.rds"))
  ## Steady-state read: relative change of the chunk means, chunk to chunk.
  rel_change <- function(x) c(NA, abs(diff(x)) / pmax(abs(head(x, -1)), 1e-300))
  stamp("relative change chunk to chunk: density %s | height_max %s | theta_1 %s",
        paste(signif(rel_change(table$density_total), 2), collapse = " "),
        paste(signif(rel_change(table$height_max), 2), collapse = " "),
        paste(signif(rel_change(table$theta_1), 2), collapse = " "))
}
stamp("done")
