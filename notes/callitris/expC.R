source("helpers.R")
PARS_OVERRIDE <- c(HEIGHT, TOLERANT)
out <- list()
## --- Thinning (Ross s3.7): 30,000/ha cohort, thinned at t = 5 (31,853/ha, ~4 m), recruits from t = 5.
T <- 25.5; N5 <- 31853
for (target in c(NA, 5000, 1000, 280, 20)) {
  ev <- if (is.na(target)) NULL else harvest(time = 5.01, fraction = 1 - target / N5, species = 1)
  r <- run_case(list(pulse_br(3, 0, 1, T), pulse_br(1, 5, T, T)), until = T, events = ev)
  ss <- stand_summary(r); lt <- ground_light(r)
  g <- function(sp, t) ss %>% filter(species == sp) %>% slice(which.min(abs(time - t)))
  row <- tibble(treat = ifelse(is.na(target), "unthinned", paste0("thin to ", target)),
    over_5.25 = g(1, 5.25)$stems_ha, over_25 = g(1, 25)$stems_ha, over_h25 = g(1, 25)$h_mean,
    light0_6 = lt$light0[which.min(abs(lt$time - 6))],
    rec_10 = g(2, 10)$stems_ha, rec_15 = g(2, 15)$stems_ha, rec_25 = g(2, 25)$stems_ha, rec_h25 = g(2, 25)$h_mean)
  out$thin[[length(out$thin) + 1]] <- row; print(as.data.frame(row), digits = 3)
}
## --- Fire and grazing: 1,000/ha overstorey sown at 0, recruits seeded from 20 (as expB).
T <- 45
base <- list(pulse_br(0.1, 0, 1, T), pulse_br(1, 20, T, T))
cases <- list(
  none = NULL,
  fire30_small = harvest(time = 30, fraction = 0.95, size_max = 1.5),
  fire30_all_sizes = harvest(time = 30, fraction = 0.95),
  graze20_30 = harvest(time = seq(20.125, 29.875, by = 0.25), fraction = 0.5, size_max = 1, species = 2))
for (k in names(cases)) {
  r <- run_case(base, until = T, events = cases[[k]])
  ss <- stand_summary(r)
  g <- function(sp, t) ss %>% filter(species == sp) %>% slice(which.min(abs(time - t)))
  row <- tibble(case = k, over_29.75 = g(1, 29.75)$stems_ha, over_30.25 = g(1, 30.25)$stems_ha, over_44 = g(1, 44.75)$stems_ha,
    rec_25 = g(2, 25)$stems_ha, rec_29.75 = g(2, 29.75)$stems_ha, rec_h29.75 = g(2, 29.75)$h_mean,
    rec_30.25 = g(2, 30.25)$stems_ha, rec_35 = g(2, 35)$stems_ha, rec_44 = g(2, 44.75)$stems_ha)
  if (!is.null(r$event_log)) row$log_rows <- NROW(r$event_log)
  out$fg[[k]] <- row; print(as.data.frame(row), digits = 3)
}
saveRDS(out, "expC.rds")
