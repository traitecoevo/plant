source("soil_cascade.R")
options(digits = 6)

source("exp2_ctl.R")   # odelia OdeControl + C++ max semantics

ck_stage <- function(y, h, rain, uptake, variant, p = P) {
  b  <- list(c(), 1/5, c(3/40, 9/40), c(3/10, -9/10, 6/5),
             c(-11/54, 5/2, -70/27, 35/27),
             c(1631/55296, 175/512, 575/13824, 44275/110592, 253/4096))
  c5 <- c(37/378, 0, 250/621, 125/594, 0, 512/1771)
  c4 <- c(2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4)
  k <- vector("list", 6); k[[1]] <- f(y, rain, uptake, variant, p)
  for (s in 2:6) {
    ys <- y; for (j in seq_len(s - 1)) ys <- ys + h * b[[s]][j] * k[[j]]
    k[[s]] <- f(ys, rain, uptake, variant, p)
  }
  y5 <- y; y4 <- y
  for (s in 1:6) { y5 <- y5 + h * c5[s] * k[[s]]; y4 <- y4 + h * c4[s] * k[[s]] }
  list(y = y5, yerr = y5 - y4, dydt_out = k[[6]])
}

cat("=== Q6. With odelia's EXACT controller + plant defaults, is the overshoot rejected? ===\n")
th0 <- c(P$theta_sat, 0.05, 0.05, 0.05, 0.05)
for (h in c(1e-4, 1e-3, 1e-2, 0.08)) {
  st <- ck_stage(th0, h, 0, rep(0, 5), "baseline")
  a  <- adjust(h, 5, st$y, st$yerr, st$dydt_out)
  cat(sprintf("  h=%7.5g: theta_1 after=%12.5g  max|yerr|=%10.3g  -> %s (h_next=%9.3g)\n",
              h, st$y[2], max(abs(st$yerr)),
              if (a$shrank) "REJECT" else "ACCEPT", a$h))
}

cat("\n=== Q7. The NaN path: what does the controller do with a non-finite error? ===\n")
y <- c(0.2, 0.2, 0.2, 0.2, 0.2)
for (case in c("finite-large", "NaN", "Inf")) {
  yerr <- rep(1e-3, 5)
  if (case == "NaN") yerr[3] <- NaN
  if (case == "Inf") yerr[3] <- Inf
  a <- adjust(1e-3, 5, y, yerr, rep(0, 5))
  cat(sprintf("  yerr=%-12s -> %s, h_next = %9.3g  (%s)\n", case,
              if (a$shrank) "REJECT" else "ACCEPT", a$h,
              if (a$h > 1e-3) sprintf("step GREW %.1fx", a$h / 1e-3) else "step held/shrank"))
}
cat("\n  => An INTERIOR NaN rejects -- but only incidentally: std::max wipes it\n")
cat("     (see cxx_max in exp2_ctl.R) and the trailing finite element rejects on\n")
cat("     its own magnitude. The NaN component itself is never accounted for.\n")
cat("     Two ways that bites, both fixed by odelia PR#54:\n")
cat("       (a) NaN in the LAST element (or all): rmax stays NaN, and since\n")
cat("           NaN > 1.1 and NaN < 0.5 are both false, control falls through\n")
cat("           to the 'no shrink' branch and the step is accepted.\n")
cat("       (b) NaN wiped by finite elements whose ratios are SMALL: rmax is\n")
cat("           finite and passes, so the step is accepted carrying a NaN.\n")
cat("     (b) is what exp4 hits. Demonstrating both:\n")
for (case in list(c(1e-3, 1e-3, 1e-3, 1e-3, NaN),
                  c(NaN, 1e-12, 1e-12, 1e-12, 1e-12))) {
  a <- adjust(1e-3, 5, rep(0.2, 5), case, rep(0, 5))
  cat(sprintf("       yerr=[%s] -> %s\n",
              paste(formatC(case, width = 6, format = "g"), collapse = ","),
              if (a$shrank) "reject" else "ACCEPT (NaN committed)"))
}

cat("\n=== Q8. Full integration of the worst case with odelia's controller ===\n")
run_odelia <- function(variant, t_end = 0.5, h_max = 5, throw_on_nonfinite = TRUE) {
  y <- th0; t <- 0; h <- CTL$h_init
  ctl <- CTL; ctl$h_max <- h_max
  ns <- 0; nr <- 0; tmax <- max(y); status <- "ok"
  while (t < t_end && ns < 2e5) {
    hh <- if (t + h > t_end) t_end - t else h
    st <- ck_stage(y, hh, 0, rep(0, 5), variant)
    if (throw_on_nonfinite && any(!is.finite(st$y))) { status <- "THROW (non-finite stage)"; break }
    a <- adjust(hh, 5, st$y, st$yerr, st$dydt_out, ctl)
    if (a$shrank) {
      nr <- nr + 1
      if (a$h >= hh) { status <- "STOP: cannot achieve accuracy"; break }
      h <- a$h
    } else {
      t <- t + hh; y <- st$y; ns <- ns + 1
      tmax <- max(tmax, max(y)); h <- a$h
    }
  }
  list(ns = ns, nr = nr, tmax = tmax, status = status, y = y)
}
for (v in c("baseline", "recv")) for (hm in c(5, 0.05)) {
  r <- run_odelia(v, h_max = hm)
  cat(sprintf("  %-9s h_max=%5.3g steps=%6d rej=%5d  max theta=%12.6g  %s\n",
              v, hm, r$ns, r$nr, r$tmax, r$status))
}

cat("\n=== Q9. Does a PULSE jump need capacity capping? (applied outside the integrator) ===\n")
cat(sprintf("  layer 0: dz = %.2f m, theta_sat = %.3f, pore space from theta_0:\n", P$dz[1], P$theta_sat))
for (t0 in c(0.10, 0.214, 0.35, 0.40)) {
  cap_m <- (P$theta_sat - t0) * P$dz[1]
  cat(sprintf("    theta_0=%.3f -> free capacity = %6.2f mm; ", t0, cap_m * 1000))
  for (d_mm in c(5, 13.3, 50)) {
    dth <- (d_mm / 1000) / P$dz[1]
    cat(sprintf(" %gmm->%.3f%s ", d_mm, t0 + dth, if (t0 + dth > P$theta_sat) "!" else ""))
  }
  cat("\n")
}
cat("  ('!' = exceeds saturation, i.e. would need capping)\n")
