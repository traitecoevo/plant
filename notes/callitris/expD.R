source("helpers.R")
LOWMORT <- list(a_dG1 = function(x) 0.1 * x, d_I = function(x) 0.001)
out <- list()
for (mv in c("tolerant", "tolerant+lowmort")) {
  PARS_OVERRIDE <- if (mv == "tolerant") c(HEIGHT, TOLERANT) else c(HEIGHT, TOLERANT, LOWMORT)
  ## Dense cohort, 50 years.
  r <- run_case(list(pulse_br(3, 0, 1, 50)), until = 50)
  ss <- stand_summary(r)
  tt <- c(5, 10, 20, 30, 49.75)
  x <- ss[sapply(tt, function(t) which.min(abs(ss$time - t))), c("time", "stems_ha", "h_mean", "h_top")]
  cat("\n==", mv, ": dense cohort (30,000/ha)\n"); print(as.data.frame(x), digits = 3)
  ## Ross s3.6: 20/ha overstorey (~10 m at t = 12), recruits from 12 to 62, wet years at 14 and 19.
  T <- 62; met <- scale_rain(MET0, c("14" = 2, "19" = 2))
  for (rec in c(TRUE, FALSE)) {
    b2 <- if (rec) pulse_br(1, 12, T, T) else pulse_br(1e-10, 12, T, T)
    r <- run_case(list(pulse_br(0.002, 0, 1, T), b2), until = T, met = met)
    ss <- stand_summary(r)
    g <- function(sp, t) ss %>% filter(species == sp) %>% slice(which.min(abs(time - t)))
    cat(sprintf("%s hist rec=%s | over@12: %.1f/ha %.1f m | over@61.75: %.1f/ha h=%.1f m | recruits@61.75: %.0f/ha hmean=%.2f\n",
      mv, rec, g(1, 12)$stems_ha, g(1, 12)$h_mean, g(1, 61.75)$stems_ha, g(1, 61.75)$h_mean,
      g(2, 61.75)$stems_ha, g(2, 61.75)$h_mean))
    if (rec) { rt <- recruit_table(r, 2, 61.75, 1); print(as.data.frame(rt %>% filter(year %in% 12:22)), digits = 3) }
  }
}
