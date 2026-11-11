## Helpers for the plant#275 Callitris experiments.
## Two "species" with identical traits stand in for two cohorts of one species:
## 1 = a cohort sown in a pulse at the start, 2 = later recruits. Seed arrival is
## controlled by each species' node schedule (seed arrives only at scheduled
## introductions), so a short schedule is a seed pulse.
suppressMessages({library(plant); library(fieldtrial); library(dplyr)})

SET <- parameter_set("callitris_central_west_nsw.yml")
MET0 <- SET$met

## Parameter overrides applied to every strategy after the trade-offs, as
## functions of the set's value. HEIGHT: hmat below the maximum height, with a
## gentler reproductive-allocation curve setting the asymptote (allocation 0.95
## at 20 m, mid Ross-2008's 16-23 m site maxima). Calibration proxies.
HEIGHT <- list(hmat = function(x) 12, a_f2 = function(x) 4.4)
## TOLERANT: uncited sensitivity variant, P50s and psi_crit doubled.
TOLERANT <- list(stem_P50 = function(x) 2 * x, root_P50 = function(x) 2 * x,
                 psi_crit = function(x) 2 * x)
PARS_OVERRIDE <- HEIGHT

drivers_all <- c("SWdown", "atm_vpd", "atm_vpd_diurnal_amp", "atm_kpa",
                 "wind_speed", "ca", "rainfall", "leaf_temp")

## Rain multiplier by calendar-year index (0-based years from t = 0).
scale_rain <- function(met, mult_by_year) {
  yr <- as.integer(format(met$date, "%Y")) - as.integer(format(met$date[1], "%Y"))
  f <- mult_by_year[as.character(yr)]
  f[is.na(f)] <- 1
  met$Precip <- met$Precip * f
  met
}

make_env <- function(met = MET0, until) {
  env <- plant::Environment("TF24")
  apply_drivers(env, met, rainfall = "spline", drivers = drivers_all, until = until)
  env
}

## A birth rate (seeds m-2 yr-1) of b between `from` and `to`, and a floor
## elsewhere; history is only recorded at introductions, so every cohort is
## introduced on the same regular grid and a pulse is carried by the birth rate.
pulse_br <- function(b, from, to, until, floor = 1e-10, edge = 1e-3) {
  x <- sort(unique(c(0, max(0, from - edge), from, to, to + edge, until + 1)))
  y <- ifelse(x >= from & x <= to, b, floor)
  list(x = x, y = y)
}

## birth: list with one entry per cohort, each a number or a pulse_br() list.
make_p <- function(birth, until, dt = 0.25) {
  n <- length(birth)
  p <- SET$p
  p$max_patch_lifetime <- until
  m <- matrix(rep(c(SET$meta$traits$lma$value, SET$meta$traits$rho$value,
                    SET$meta$traits$omega$value), each = n), nrow = n,
              dimnames = list(NULL, c("lma", "rho", "omega")))
  p <- plant::add_strategies(p, m, SET$tradeoffs, birth_rate = birth,
                             keep_existing = FALSE)
  p$hyperpar <- SET$tradeoffs
  for (i in seq_along(p$strategies)) for (k in names(PARS_OVERRIDE))
    p$strategies[[i]]$pars[[k]] <- PARS_OVERRIDE[[k]](p$strategies[[i]]$pars[[k]])
  p$max_patch_lifetime <- until
  sched <- seq(0, until - dt / 2, by = dt)
  p$node_schedule_times <- rep(list(sched), n)
  p
}

run_case <- function(birth, until, met = MET0, events = NULL,
                     ctrl = control(), dt = 0.25) {
  p <- make_p(birth, until, dt)
  env <- make_env(met, until)
  ev <- if (is.null(events)) NULL else plant::events(plant::events_default(p), events)
  t0 <- Sys.time()
  res <- run_scm(p, env, ctrl, collect = TRUE, events = ev)
  attr(res, "elapsed") <- as.numeric(Sys.time() - t0, units = "secs")
  res
}

## Per time and cohort: stems ha-1, top and density-weighted mean height, and
## basal area (m2 ha-1) from stem area if available.
stand_summary <- function(res) {
  sp <- res$species
  dens <- integrate_over_size_distribution(sp) %>%
    transmute(time, species, stems_ha = 1e4 * density,
              h_mean = height / density, storage_mean = storage / density,
              npp_mean = net_mass_production_dt / density)
  ## Top height: tallest node still carrying a non-negligible share of stems.
  hts <- sp %>% filter(!is.na(height), step > 1) %>%
    group_by(time, species) %>%
    summarise(h_top = max(height), .groups = "drop")
  left_join(dens, hts, by = c("time", "species"))
}

## Values at chosen times (nearest step).
at_times <- function(df, times) {
  df %>% group_by(species) %>%
    reframe(across(everything(), ~ .x[sapply(times, function(t) which.min(abs(time - t)))])) %>%
    ungroup()
}

## Recruits by birth date: density at birth per seed (an establishment index),
## and stems ha-1 surviving at time t_obs, summed per birth year.
add_birth <- function(sp) sp %>% filter(!is.na(height)) %>% group_by(species, node) %>%
  mutate(birth_time = min(time)) %>% ungroup()
recruit_table <- function(res, sp_i, t_obs, br) {
  sp <- add_birth(res$species) %>% filter(species == sp_i)
  births <- sp %>% filter(time == birth_time) %>% transmute(birth_time, n_birth = exp(log_density_state))
  s <- sp %>% filter(abs(time - t_obs) == min(abs(time - t_obs))) %>% arrange(birth_time)
  ## trapezium per segment, attributed to the segment's start year
  seg <- tibble(b0 = head(s$birth_time, -1), b1 = tail(s$birth_time, -1),
                d = 0.5 * (head(exp(s$log_density_state), -1) + tail(exp(s$log_density_state), -1)) * diff(s$birth_time))
  surv <- seg %>% group_by(year = floor(b0)) %>% summarise(stems_ha = 1e4 * sum(d))
  est <- births %>% mutate(est = n_birth / br) %>% group_by(year = floor(birth_time)) %>%
    summarise(est_index = mean(est))
  full_join(est, surv, by = "year")
}
ground_light <- function(res) res$env$light_availability %>% group_by(step, time) %>%
  slice_min(height, n = 1, with_ties = FALSE) %>% ungroup() %>% select(time, light0 = light_availability)
