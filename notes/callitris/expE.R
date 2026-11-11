source("helpers.R")
PARS_OVERRIDE <- c(HEIGHT, TOLERANT)
T <- 35; WET <- c(22, 27)
set.seed(275)
## Intermittent rain: Poisson storms, 40 a year averaging 14.3 mm (572 mm/yr); wet years double the storm rate.
pl <- do.call(rbind, lapply(0:(T - 1), function(y) {
  x <- fieldtrial::pulse_series_poisson(y, y + 1, rate = if (y %in% WET) 80 else 40, mean_amount = 0.0143)
  x }))
cat("annual rain (mm):", round(tapply(pl$amount, floor(pl$time), sum) * 1000)[c(1, 20:29)], "\n")
for (n0 in c(0, 1000)) {
  p <- make_p(list(pulse_br(max(n0, 1e-6) / 1e4, 0, 1, T), pulse_br(1, 20, T, T)), T)
  env <- plant::Environment("TF24")
  apply_drivers(env, MET0, rainfall = "spline", drivers = setdiff(drivers_all, "rainfall"), until = T)
  env$extrinsic_drivers_set_constant("rainfall", 0)
  ev <- plant::events(plant::events_default(p), plant::rainfall_pulse(pl$time, pl$amount, env = env))
  r <- run_scm(p, env, control(), collect = TRUE, events = ev)
  rt <- recruit_table(r, 2, T - 0.25, 1) %>% filter(year >= 20)
  e <- r$env$soil_moist %>% group_by(step) %>% mutate(layer = row_number()) %>% filter(layer == 1) %>% ungroup()
  e$psi <- sapply(e$soil_moist, env$psi_from_soil_moist)
  ps <- e %>% filter(time >= 20) %>% group_by(year = floor(time)) %>% summarise(psi_top_median = median(psi), psi_top_max = max(psi))
  cat("\n== Poisson rain, overstorey", n0, "/ha\n"); print(as.data.frame(left_join(rt, ps, by = "year")), digits = 3)
}
