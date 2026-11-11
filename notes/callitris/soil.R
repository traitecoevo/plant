source("helpers.R")
for (case in list(list("default", HEIGHT, 20), list("default", HEIGHT, 1000), list("tolerant", c(HEIGHT, TOLERANT), 30000))) {
  PARS_OVERRIDE <- case[[2]]
  p <- make_p(list(pulse_br(case[[3]] / 1e4, 0, 1, 30)), 30); env <- make_env(MET0, 30)
  r <- run_scm(p, env, control(), collect = TRUE)
  e <- r$env$soil_moist %>% group_by(step) %>% mutate(layer = row_number()) %>% ungroup()
  sat <- env$soil_moist_sat
  cat("\n==", case[[1]], case[[3]], "stems/ha; soil_moist_sat", sat, "; widths", round(env$get_soil_layer_widths(), 3), "\n")
  s <- e %>% mutate(period = cut(time, c(-1, 5, 15, 30))) %>% group_by(period, layer) %>%
    summarise(mean = mean(soil_moist), min = min(soil_moist), .groups = "drop")
  s$psi_min <- mapply(function(m, l) env$psi_from_soil_moist(m), s$min, s$layer)
  print(as.data.frame(s), digits = 3)
}
