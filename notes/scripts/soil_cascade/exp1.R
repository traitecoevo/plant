source("soil_cascade.R")
options(digits = 6)

cat("=== Q1. Is theta_sat a barrier in the EXACT vector field? ===\n")
cat("Worst case: layer above saturated, layer below dry.\n")
for (tb in c(0.05, 0.2, 0.35, 0.9 * P$theta_sat, 0.99 * P$theta_sat, P$theta_sat)) {
  th <- c(P$theta_sat, tb, rep(0.05, 3))
  r <- rates(th, rain = 0, uptake = rep(0, 5), "baseline")
  cat(sprintf("  theta_1=%.5f (%.3f of sat)  in=%8.3f out=%8.3f  rate_1=%10.3f /yr\n",
              tb, tb / P$theta_sat, r$inflow[2], K_of(tb), r$rate[2]))
}
cat("\n  At theta = theta_sat exactly, both in and out equal K_sat, so rate = 0.\n")
cat("  => the field cannot push theta above saturation. Overshoot must be NUMERICAL.\n")

cat("\n=== Q2. Reproduce the blow-up: explicit step across the boundary layer ===\n")
th0 <- c(P$theta_sat, 0.05, 0.05, 0.05, 0.05)
for (h in c(1e-5, 1e-4, 1e-3, 1e-2, 0.08)) {
  r <- f(th0, 0, rep(0, 5), "baseline")
  th1 <- th0 + h * r
  cat(sprintf("  h=%7.5g yr:  rate_1=%8.2f  ->  theta_1=%12.5f  (%s)\n",
              h, r[2], th1[2],
              if (th1[2] > P$theta_sat) sprintf("PAST SAT by %.3g", th1[2] - P$theta_sat) else "ok"))
}
cat(sprintf("\n  Stability limit for an explicit step: h < 2/|lambda_max|.\n"))

cat("\n=== Q3. Jacobian structure and spectrum ===\n")
probe <- c(0.42, 0.40, 0.30, 0.20, 0.10)
for (v in c("baseline", "recv", "donor")) {
  J <- jacobian(probe, rain = 1, uptake = rep(0, 5), v)
  upper <- max(abs(J[upper.tri(J)]))
  ev <- eigen(J, only.values = TRUE)$values
  cat(sprintf("  %-9s max|upper-tri| = %9.3g   max|lambda| = %10.4g   lower-tri: %s\n",
              v, upper, max(abs(ev)), if (upper < 1e-6) "YES" else "NO"))
}

cat("\n  Spectrum near saturation (theta -> theta_sat), baseline vs recv:\n")
for (frac in c(0.5, 0.9, 0.99, 1.0)) {
  th <- rep(frac * P$theta_sat, 5)
  eb <- max(abs(eigen(jacobian(th, 1, rep(0, 5), "baseline"), only.values = TRUE)$values))
  er <- max(abs(eigen(jacobian(th, 1, rep(0, 5), "recv"), only.values = TRUE)$values))
  cat(sprintf("    theta/sat=%.2f   baseline |lambda|=%10.4g   recv |lambda|=%10.4g   ratio %.2fx\n",
              frac, eb, er, er / eb))
}

cat("\n=== Q4. Does the receiver limiter reduce the arrival rate? ===\n")
cat("  (inflow into a partly-wet layer beneath a saturated one)\n")
for (frac in c(0.5, 0.9, 0.99)) {
  th <- c(P$theta_sat, frac * P$theta_sat, rep(0.05, 3))
  rb <- rates(th, 0, rep(0, 5), "baseline")
  rr <- rates(th, 0, rep(0, 5), "recv")
  cat(sprintf("    theta_1/sat=%.2f  baseline rate=%9.3f  recv rate=%9.3f  (%.0f%% of baseline)\n",
              frac, rb$rate[2], rr$rate[2], 100 * rr$rate[2] / rb$rate[2]))
}

cat("\n=== Q5. Integrate the worst case: does it survive, and at what cost? ===\n")
for (v in c("baseline", "recv")) {
  for (hmax in c(5, 0.05)) {
    r <- ck(th0, t_end = 0.5, rain = 0, uptake = rep(0, 5), variant = v, h_max = hmax)
    cat(sprintf("  %-9s h_max=%5.3g  steps=%7d rej=%6d  max theta seen=%12.5g  finite=%s\n",
                v, hmax, r$nstep, r$nrej, r$theta_max, r$finite))
  }
}
