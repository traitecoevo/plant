## Shared by fixed_point.R and long_patch.R (#519): settings from the
## environment, logging, and the pulsed TF24 stand the long-patch tasks share.
## Sourced, not run. Every setting is an environment variable named PLANT_519_*,
## so one script serves a smoke test and a full run; run from the plant root
## with the library that holds the pinned odelia and phylloptim:
##   R_LIBS=<lib> Rscript --no-init-file notes/scripts/continuous-patch-519/<script>.R

suppressPackageStartupMessages(library(plant))

num_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(v)
}
vec_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(strsplit(v, ",")[[1]])
}
stamp <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep = " ")
  flush(stdout())
}
## Which of a script's tasks to run: PLANT_519_TASK, comma-separated, or the
## script's default.
tasks_to_run <- function(all, default = all) {
  v <- Sys.getenv("PLANT_519_TASK", "")
  chosen <- if (v == "") default else strsplit(v, ",")[[1]]
  unknown <- setdiff(chosen, all)
  if (length(unknown) > 0) {
    stop("unknown task ", paste(unknown, collapse = ", "), "; tasks are ",
         paste(all, collapse = ", "))
  }
  chosen
}
out_dir <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
trapezium <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

## ---- The pulsed TF24 stand ------------------------------------------------
## One TF24 patch, no disturbance, rainfall as a stationary marked Poisson
## series of pulses (PLANT_519_LAMBDA a year, exponential depths of mean
## PLANT_519_DEPTH m), the continuous driver off, a cohort a year, and the top
## layer at PLANT_519_ABSORB. A seed fixes the series: one realisation shared by
## every run compared on it.
pulsed <- list(
  lambda      = num_env("PLANT_519_LAMBDA", 20),
  depth       = num_env("PLANT_519_DEPTH", 0.014),
  series_tmax = num_env("PLANT_519_SERIES_TMAX", 1000),
  absorb      = num_env("PLANT_519_ABSORB", 0.05),
  resident    = 0.0825)

pulse_series_for <- function(seed) {
  set.seed(seed)
  pulse_series_poisson(0, pulsed$series_tmax, rate = pulsed$lambda,
                       mean_amount = pulsed$depth)
}
## The series' pulses in (a, b], as events.
pulses_between <- function(series, a, b) {
  k <- series$time > a & series$time <= b
  if (any(k)) rainfall_pulse(series$time[k], series$amount[k])
}
pulsed_parameters <- function(lma, lifetime, birth = 1) {
  p <- scm_base_parameters("TF24")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- lifetime
  add_strategies(p, trait_matrix(lma, "lma"), hyperpar = TF24_hyperpar,
                 birth_rate = birth)
}
## TF24's soil ODE is stiff for an explicit stepper without a cap on the step.
pulsed_control <- function(absorb = pulsed$absorb) {
  control(node_density_in_birth_date = TRUE, ode_step_size_max = 0.05,
          absorb_epsilon = absorb)
}
pulsed_env <- function() {
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 0)
  env
}
