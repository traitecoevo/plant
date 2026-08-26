# Standalone replica of TF24's soil-water cascade, for the Stage A design pass.
# Mirrors tf24_environment.h compute_rates() / soil_K_from_soil_theta() exactly.

P <- list(
  n         = 5L,
  depth     = 1.5,
  theta_sat = 0.428,
  K_sat     = 163.0411,
  n_psi     = 6.57,
  a_infil   = 1,
  b_infil   = 8,
  theta_res = 1e-2
)
P$dz  <- rep(P$depth / P$n, P$n)
P$b_K <- 2 * P$n_psi + 3        # 16.14

# tf24_environment.h:403 -- theta clamped to [0, theta_sat] inside K
K_of <- function(theta, p = P) {
  t <- pmin(pmax(theta, 0), p$theta_sat)
  p$K_sat * (t / p$theta_sat)^p$b_K
}

# The saturation-excess factor used at the surface (line 336).
excess_factor <- function(theta, p = P) {
  pmax(0, 1 - p$a_infil * (theta / p$theta_sat)^p$b_infil)
}

# variant: "baseline" = donor-only (as shipped)
#          "recv"     = receiver rejects excess, water leaves as rejected flux
#          "donor"    = donor throttles when receiver is full (breaks triangularity)
rates <- function(theta, rain, uptake, variant = "baseline", p = P) {
  n <- p$n
  out <- K_of(theta, p)                      # gravity drainage, donor-only
  inflow <- numeric(n)
  rejected <- numeric(n)

  inflow[1] <- rain * excess_factor(theta[1], p)

  for (i in seq_len(n)) {
    if (i > 1) {
      avail <- out[i - 1]
      if (variant == "baseline") {
        inflow[i] <- avail
      } else if (variant == "recv") {
        inflow[i]   <- avail * excess_factor(theta[i], p)
        rejected[i] <- avail - inflow[i]
      } else if (variant == "donor") {
        inflow[i] <- avail * excess_factor(theta[i], p)
      }
    }
  }
  # donor variant: the donor's own outflow is throttled by the receiver, so the
  # donor's rate must use the throttled value -> rate_{i-1} depends on theta_i.
  out_eff <- out
  if (variant == "donor") {
    for (i in 2:n) out_eff[i - 1] <- inflow[i]
  }

  r <- (inflow - out_eff - uptake) / p$dz

  # positivity guard, tf24_environment.h:375-377
  guard <- theta <= p$theta_res & !(r > 0)
  r[guard] <- 0
  list(rate = r, rejected = rejected, drain = out_eff[n], inflow = inflow)
}

f <- function(theta, rain, uptake, variant, p = P) {
  rates(theta, rain, uptake, variant, p)$rate
}

jacobian <- function(theta, rain, uptake, variant, p = P, h = 1e-8) {
  n <- length(theta)
  J <- matrix(0, n, n)
  f0 <- f(theta, rain, uptake, variant, p)
  for (j in seq_len(n)) {
    tp <- theta; tp[j] <- tp[j] + h
    J[, j] <- (f(tp, rain, uptake, variant, p) - f0) / h
  }
  J
}

# Cash-Karp adaptive RK, matching odelia's controller shape closely enough to
# count steps and detect domain exit.
ck <- function(theta0, t_end, rain, uptake, variant, p = P,
               h_max = 5, h_init = 1e-3, tol = 1e-6) {
  a  <- c(0, 1/5, 3/10, 3/5, 1, 7/8)
  b  <- list(c(), 1/5, c(3/40, 9/40), c(3/10, -9/10, 6/5),
             c(-11/54, 5/2, -70/27, 35/27),
             c(1631/55296, 175/512, 575/13824, 44275/110592, 253/4096))
  c5 <- c(37/378, 0, 250/621, 125/594, 0, 512/1771)
  c4 <- c(2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4)

  y <- theta0; t <- 0; h <- h_init
  nstep <- 0; nrej <- 0; theta_max <- max(y); hmin_seen <- Inf
  while (t < t_end && nstep < 2e6) {
    if (t + h > t_end) h <- t_end - t
    k <- vector("list", 6)
    k[[1]] <- f(y, rain, uptake, variant, p)
    for (s in 2:6) {
      ys <- y
      for (j in seq_len(s - 1)) ys <- ys + h * b[[s]][j] * k[[j]]
      k[[s]] <- f(ys, rain, uptake, variant, p)
    }
    y5 <- y; y4 <- y
    for (s in 1:6) { y5 <- y5 + h * c5[s] * k[[s]]; y4 <- y4 + h * c4[s] * k[[s]] }
    err <- max(abs(y5 - y4) / (abs(y) + 1e-10)) / tol
    if (!is.finite(err) || err > 1.1) {
      nrej <- nrej + 1
      h <- max(h * max(0.2, 0.9 / err^0.2), 1e-14)
    } else {
      t <- t + h; y <- y5; nstep <- nstep + 1
      theta_max <- max(theta_max, max(y))
      hmin_seen <- min(hmin_seen, h)
      h <- min(h * min(5, 0.9 / max(err, 1e-10)^0.2), h_max)
    }
    if (any(!is.finite(y))) break
  }
  list(theta = y, t = t, nstep = nstep, nrej = nrej,
       theta_max = theta_max, h_min = hmin_seen, finite = all(is.finite(y)))
}
