## What a window run costs, and which cohorts it integrates (#519): the pulsed
## TF24 stand of landscape.R at the resident's equilibrium birth rate on series 1,
## burnt in to 200 years, exported, and resumed into the window [200, 300] run to
## 350, timed in 25-year legs with the cohorts still integrated at each.
##
##   PLANT_519_MODE    "absorb" (this branch: the top layer at epsilon 0.05) or
##                     "floor", which needs a build of c706e257, the commit with
##                     the survival floor at 1e-9 that the top layer replaced;
##                     on any later build it errors. For that build the script
##                     builds its SCMs with plant:::SCM() and falls back to
##                     window_regime(), which it predates run_scm(run = FALSE)
##                     and fitness_window() by.
##   PLANT_519_BIRTH   resident birth rate (8.17503, b* on series 1)
##   PLANT_519_OUT     output directory

suppressPackageStartupMessages(library(plant))
mode  <- Sys.getenv("PLANT_519_MODE", "absorb")
birth <- as.numeric(Sys.getenv("PLANT_519_BIRTH", "8.17503"))
out_dir <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(format(Sys.time(), "%H:%M:%S"), "plant", as.character(packageVersion("plant")), "at",
    find.package("plant"), "; mode", mode, "; birth rate", birth, "\n")

## The series and stand of workflow.R and landscape.R, seed 1.
set.seed(1)
series_tmax <- 1000
## pulse_series_poisson() draws count, times, then amounts, so a seed rebuilds
## the series this script's committed numbers were measured on. The floor-era
## build predates it, so the same draws are spelled out for that one.
draw_pulses <- if (exists("pulse_series_poisson", mode = "function")) pulse_series_poisson else
  function(from, to, rate, mean_amount) {
    n <- rpois(1, rate * (to - from))
    data.frame(time = sort(runif(n, from, to)), amount = rexp(n, rate = 1 / mean_amount))
  }
series      <- draw_pulses(0, series_tmax, rate = 20, mean_amount = 0.014)
pulse_time  <- series$time
pulse_depth <- series$amount
n_pulses    <- nrow(series)
intro_times <- seq(0, series_tmax, by = 1)
ctrl <- if (mode == "floor") {
  control(node_density_in_birth_date = TRUE, ode_step_size_max = 0.05, retire_survival = 1e-9)
} else {
  control(node_density_in_birth_date = TRUE, ode_step_size_max = 0.05, absorb_epsilon = 0.05)
}
base_parameters <- function(lifetime) {
  p <- scm_base_parameters("TF24")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- lifetime
  add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = TF24_hyperpar, birth_rate = birth)
}
make_env <- function() {
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 0)
  env
}
pulses_in <- function(a, b) {
  k <- pulse_time > a & pulse_time <= b
  rainfall_pulse(pulse_time[k], pulse_depth[k])
}

p <- base_parameters(200)
p$node_schedule_times <- list(intro_times[intro_times < 200])
scm <- plant:::SCM("TF24", "TF24_Env")(p, make_env(), events(node_introductions(p), pulses_in(0, 200)), ctrl)
w_burn <- system.time(scm$run_to(200))[["elapsed"]]
state <- export_patch_state(scm)

p2 <- base_parameters(350)
p2 <- set_initial_state(p2, state)
## The floor-era build names it window_regime().
weight_window <- if (exists("fitness_window", mode = "function")) fitness_window else get("window_regime")
p2 <- weight_window(p2, 200, 300)
p2$node_schedule_times <- list(intro_times[intro_times > max(state$node_times[[1]]) & intro_times < 350])
scm2 <- plant:::SCM("TF24", "TF24_Env")(p2, make_env(), events(node_introductions(p2), pulses_in(200, 350)), ctrl)
rows <- lapply(seq(225, 350, by = 25), function(t) {
  w <- system.time(scm2$run_to(t))[["elapsed"]]
  sp <- scm2$patch$species[[1]]
  in_window <- sp$node_times >= 200 & sp$node_times <= 300
  data.frame(mode = mode, time = t, wall_s = w, steps = length(scm2$ode_times), nodes = sp$size,
             integrated = sp$size - sum(sp$retired), window_integrated = sum(in_window & !sp$retired))
})
d <- do.call(rbind, rows)
d$burn_in_s <- w_burn
d$R_window <- scm2$net_reproduction_ratios
print(d, row.names = FALSE)
cat(sprintf("%s: burn-in %.1f s, window run %.1f s, %d steps, R %.6f\n", mode, w_burn, sum(d$wall_s),
            max(d$steps), scm2$net_reproduction_ratios))
write.csv(d, file.path(out_dir, paste0("window-cost-", mode, ".csv")), row.names = FALSE)
