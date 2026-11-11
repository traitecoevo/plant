source("helpers.R")
diag_run <- function(label, met = MET0, until = 40, set_leaf25 = FALSE, p_mod = NULL) {
  p <- make_p(list(pulse_br(0.002, 0, 1, until)), until)
  if (!is.null(p_mod)) p <- p_mod(p)
  env <- plant::Environment("TF24")
  drv <- if (set_leaf25) setdiff(drivers_all, "leaf_temp") else drivers_all
  apply_drivers(env, met, rainfall = "spline", drivers = drv, until = until)
  r <- run_scm(p, env, control(), collect = TRUE)
  ss <- stand_summary(r)
  sp <- r$species %>% filter(!is.na(height), step > 1) %>% group_by(time) %>%
    slice_max(height, n = 1) %>% ungroup() %>%
    select(time, height, assimilation, transpiration, opt_psi_leaf, storage, net_mass_production_dt)
  en <- r$env
  cat("\n==", label, "\n")
  x <- at_times(ss, c(5, 10, 20, 30, 39.5)) %>% select(time, stems_ha, h_top)
  y <- sp[sapply(c(5, 10, 20, 30, 39.5), function(t) which.min(abs(sp$time - t))), ]
  print(as.data.frame(cbind(x, y[, -1:-2])), digits = 3)
  invisible(r)
}
r <- diag_run("baseline (572 mm, leaf_temp = noon Tair)")
print(names(r$env)); 
diag_run("leaf_temp fixed at plant default 25 C", set_leaf25 = TRUE)
diag_run("rain x3", met = scale_rain(MET0, setNames(rep(3, 106), 0:105)))
diag_run("P50s x3 (stem, root) -- sensitivity only", p_mod = function(p) {
  for (i in seq_along(p$strategies)) { p$strategies[[i]]$pars$stem_P50 <- 3 * p$strategies[[i]]$pars$stem_P50;
     p$strategies[[i]]$pars$root_P50 <- 3 * p$strategies[[i]]$pars$root_P50;
     p$strategies[[i]]$pars$psi_crit <- 3 * p$strategies[[i]]$pars$psi_crit }; p })
