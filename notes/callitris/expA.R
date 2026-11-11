source("helpers.R")
T <- 60; out <- list()
for (hyd in c("default", "tolerant")) {
  PARS_OVERRIDE <- if (hyd == "default") HEIGHT else c(HEIGHT, TOLERANT)
  for (n_ha in c(20, 1000, 10000, 30000)) {
    r <- run_case(list(pulse_br(n_ha / 1e4, 0, 1, T)), until = T)
    out[[paste(hyd, n_ha)]] <- stand_summary(r) %>% mutate(hyd = hyd, n0 = n_ha, rain = 1)
  }
  for (f in c(0.8, 1.25)) {
    met <- scale_rain(MET0, setNames(rep(f, 106), 0:105))
    r <- run_case(list(pulse_br(0.1, 0, 1, T)), until = T, met = met)
    out[[paste(hyd, "rain", f)]] <- stand_summary(r) %>% mutate(hyd = hyd, n0 = 1000, rain = f)
  }
}
A <- bind_rows(out); saveRDS(A, "expA.rds")
tt <- c(5, 10, 20, 30, 45, 59.5)
A %>% group_by(hyd, n0, rain) %>%
  reframe(t = tt, stems = stems_ha[sapply(tt, function(x) which.min(abs(time - x)))],
          h_mean = h_mean[sapply(tt, function(x) which.min(abs(time - x)))],
          h_top = h_top[sapply(tt, function(x) which.min(abs(time - x)))]) %>%
  mutate(across(c(stems, h_mean, h_top), ~ signif(.x, 3))) %>%
  tidyr::pivot_wider(names_from = t, values_from = c(stems, h_mean)) %>% select(-h_top) %>%
  as.data.frame() %>% print(width = 250)
