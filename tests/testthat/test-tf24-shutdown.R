# R-C / #55: a shut-down TF24 leaf (stem held at psi_crit, transpiration impossible)
# draws no water. Before the fix, set_shutdown_state left soil_consumption_/E_up_
# stale, so a reused leaf fed a soil-moisture-independent phantom uptake into the
# soil balance and zeroed the reverse-mode gradient across the drought regime.

mkleaf <- function() Leaf(vcmax_25 = 100, jmax_25 = 100 * 167, c = 2.04, b = 3, psi_crit = 5,
  root_c = 2.65, root_b = 1.29, root_psi_crit = 1.29 * (log(1 / 0.05))^(1 / 2.65), beta2 = 1,
  a = 0.3, curv_fact_elec_trans = 0.7, curv_fact_colim = 0.99, GSS_tol_abs = 1e-8,
  vulnerability_curve_ncontrol = 100, ci_abs_tol = 1e-8, ci_niter = 1000, g1_TF24 = 46.32995,
  beta_R_H = 3.4e3, beta_R_V = 9.4e4)
setp <- function(L, psi) L$set_physiology(area_leaf = 1, mass_root_prop = rep(1 / length(psi), length(psi)),
  rho = 608, a_bio = 0.0245, PPFD = 1800, psi_soil = psi, soil_depth = seq(0.3, by = 0.3, length.out = length(psi)),
  leaf_specific_conductance_max = 1e-4, atm_vpd = 1, ca = 40, sapwood_volume_per_leaf_area = 1e-4,
  leaf_temp = 25, atm_o2_kpa = 21, atm_kpa = 101.3)

test_that("a shut-down TF24 leaf draws no water (#55)", {
  L <- mkleaf()
  setp(L, rep(6, 5))          # every layer drier than psi_crit (= 5) -> shutdown
  L$find_root_collar_psi()
  expect_true(is.finite(L$profit_))
  expect_identical(L$E_up_, 0)
  expect_identical(L$soil_consumption_, rep(0, 5))
})

test_that("a reused leaf does not carry stale uptake into shutdown (#55)", {
  L <- mkleaf()
  setp(L, rep(1, 5)); L$find_root_collar_psi()          # responsive solve leaves nonzero soil_consumption_
  expect_gt(sum(L$soil_consumption_), 0)
  setp(L, rep(6, 5)); L$find_root_collar_psi()          # shutdown on the SAME leaf must reset to 0
  expect_identical(L$E_up_, 0)
  expect_identical(L$soil_consumption_, rep(0, 5))
})

test_that("TF24 stand uptake vanishes in deep drought without a phantom floor (#55)", {
  SAT <- 0.428; KSAT <- 163.0411; DZmm <- 300; P <- 2 * 6.57 + 3
  Kf <- function(th) KSAT * (pmax(th, 1e-9) / SAT)^P
  p0 <- scm_base_parameters("TF24"); p0$max_patch_lifetime <- 30
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"))
  env <- Environment("TF24"); env$extrinsic_drivers_set_variable("rainfall", 0:2, c(0, 0, 0))
  h <- c(14, 10, 7, 5, 3.5, 2.2, 1.3, 0.7); d <- c(0.015, 0.03, 0.06, 0.12, 0.25, 0.6, 1.5, 4.0)
  st <- make_initial_state(p1, heights = h, densities = d, env = env, ctrl = Control())
  patch <- SCM("TF24", "TF24_Env")(set_initial_state(p1, st), env, Control())$patch
  patch$compute_environment()
  ne <- patch$environment$ode_size; ns <- patch$ode_size; si <- (ns - ne + 1):(ns - ne + 5)
  # aggregate per-layer uptake backed out of the real patch derivs at uniform theta (zero inflow)
  agg_uptake <- function(theta) {
    y <- patch$ode_state; y[si] <- pmin(pmax(rep(theta, 5), 1e-6), SAT - 1e-6)
    r <- patch$derivs(y, 0)[si]; K <- Kf(theta)
    up <- numeric(5); up[1] <- -K - DZmm / 1000 * r[1]  # dz in m for this backing-out is immaterial to sign/zero
    for (i in 2:5) up[i] <- -DZmm / 1000 * r[i]
    up
  }
  expect_gt(sum(abs(agg_uptake(0.20))), 1e-4)             # responsive when wet
  expect_lt(max(abs(agg_uptake(0.04))), 1e-9)            # ~0 in deep drought -- no phantom floor
  # the drought-response gradient is live in the responsive regime (was dead before the fix)
  slope <- (sum(agg_uptake(0.17)) - sum(agg_uptake(0.15))) / 0.02
  expect_gt(abs(slope), 1e-3)
})
