source("helpers.R")
PARS_OVERRIDE <- c(HEIGHT, TOLERANT)
T <- 45; WET <- c(22, 27); BR2 <- 1
wet_met <- scale_rain(MET0, setNames(rep(2, length(WET)), WET))
res <- list()
for (n0 in c(0, 1000, 30000)) for (rain in c("normal", "wet")) for (rec in c(TRUE, FALSE)) {
  if (!rec && rain == "wet") next
  b1 <- pulse_br(max(n0, 1e-6) / 1e4, 0, 1, T)
  b2 <- if (rec) pulse_br(BR2, 20, T, T) else pulse_br(1e-10, 20, T, T)
  r <- run_case(list(b1, b2), until = T, met = if (rain == "wet") wet_met else MET0)
  key <- paste(n0, rain, rec)
  res[[key]] <- list(n0 = n0, rain = rain, rec = rec, ss = stand_summary(r),
                     rt = if (rec) recruit_table(r, 2, T - 0.25, BR2) else NULL,
                     light = ground_light(r), t = attr(r, "elapsed"))
  cat(key, "elapsed", round(attr(r, "elapsed"), 1), "\n")
}
saveRDS(res, "expB.rds")
