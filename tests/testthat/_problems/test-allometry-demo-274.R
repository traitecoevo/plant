# Extracted from test-allometry-demo.R:274

# prequel ----------------------------------------------------------------------
demo_helpers <- test_path("..", "..", "overstorey_staging",
                          "allometry_demo_helpers.R")
tf24_departure_rates_r <- function(s, height, phi, psi, tracked, P, A) {
  p <- s$pars
  storage_gate_width    <- 0.1    # TF24_Strategy::storage_gate_width
  storage_prod_eps      <- 1e-4   # TF24_Strategy::storage_prod_eps
  plasticity_gate_width <- 0.02   # TF24_Strategy::plasticity_gate_width
  eta_c <- 1 - 2 / (1 + p$eta) + 1 / (1 + 2 * p$eta)   # CanopyShape::eta_c

  ## The gate reads the TRACKED marginal balance, and its effect tapers to
  ## nothing at the canopy floor.
  rho_carbon <- 1 - p$a_pl0 / (1 + exp((tracked - p$a_pl1) / plasticity_gate_width))
  shed_gate  <- 1 - exp((log(p$a_pl3) - phi) / p$a_pl2)
  rho_l <- 1 - (1 - rho_carbon) * shed_gate
  rho_s <- rho_l * exp(-psi / p$a_pl2)
  u_l   <- 1 - rho_l
  u_s   <- 1 - rho_s

  ## Reserves still gate growth, as they did before.
  A_s   <- p$theta * A * exp(psi)
  S_max <- p$a_st1 * A_s * height * eta_c * p$rho
  ## (r is supplied by the caller through P's own individual; recomputed here
  ## only for the growth gate.)
  list(rho_l = rho_l, u_l = u_l, u_s = u_s, S_max = S_max,
       finish = function(r) {
         G     <- 1 / (1 + exp(-(r - p$a_st2) / storage_gate_width))
         Ppos  <- 0.5 * (P + sqrt(P^2 + storage_prod_eps^2))
         F     <- Ppos * G
         f_g   <- 1 - p$a_f1 / (1 + exp(p$a_f2 * (1 - height / p$hmat)))
         sigma <- rho_l * (1 - exp(phi / p$a_pl2))
         c_reb <- p$lma + p$a_r1 +
           p$theta * height * eta_c * p$rho * (1 + p$a_b1)
         b <- sigma * F * f_g / (c_reb * A)
         list(dphi = b - u_l * p$k_l,
              dpsi = u_l * p$k_l - u_s * p$k_s + b * (exp(-psi) - 1))
       })
}

# test -------------------------------------------------------------------------
skip_if_not(file.exists(demo_helpers), "demo helpers not present")
skip_on_cran()
source(demo_helpers, local = TRUE)
ctl <- lapply(c(0.30, 0.20, 0.15), allometry_sapwood_controller)
names(ctl) <- c("wet", "moist", "dry")
for (nm in c("wet", "moist")) {
    d <- ctl[[nm]]
    expect_equal(d$huber_ratio[which.min(abs(d$R_s))],
                 d$huber_ratio[which.max(d$dheight)])
  }
peak <- vapply(ctl, function(d) d$huber_ratio[which.max(d$dheight)], numeric(1))
expect_equal(unname(round(peak[["wet"]], 2)), 1.28)
expect_equal(unname(round(peak[["moist"]], 2)), 1.35)
expect_gt(peak[["dry"]], peak[["moist"]])
expect_lt(ctl$dry$dheight[[1]], 1e-8)
