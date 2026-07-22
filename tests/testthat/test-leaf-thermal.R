# Leaf thermal-damage / acclimation layer (ATLS, Lumry-Eyring; issue #566).
# Phase 1: the Leaf-level damage factor N, its wiring into jmax, the constitutive
# tolerance (T_opt via d_S, T_crit_0) and repair (k_r1_0) axes, and the
# numerical-stability contract (finite, bounded, smooth, sane on outlandish input).

mk_leaf <- function() {
  theta <- 0.000157; h <- 5; K_s <- 1
  root_c <- 2.65; root_b <- 1.29
  Leaf(vcmax_25 = 100, jmax_25 = 100 * 1.67, c = 2.04, b = 3, psi_crit = 5,
       root_c = root_c, root_b = root_b,
       root_psi_crit = root_b * (log(1 / 0.05))^(1 / root_c),
       beta2 = 1, a = 0.3, curv_fact_elec_trans = 0.7, curv_fact_colim = 0.99,
       GSS_tol_abs = 1e-8, vulnerability_curve_ncontrol = 100, ci_abs_tol = 1e-6,
       ci_niter = 1000, g1_TF24 = 46.32995, beta_R_H = 3.4e3, beta_R_V = 9.4e4)
}

set_phys <- function(l, Tleaf) {
  theta <- 0.000157; h <- 5; K_s <- 1
  l$set_physiology(area_leaf = 0.05, mass_root_prop = 1, rho = 608, a_bio = 0.0245,
    PPFD = 900, psi_soil = 2, soil_depth = 1,
    leaf_specific_conductance_max = K_s * theta / h, atm_vpd = 2, ca = 40,
    sapwood_volume_per_leaf_area = theta * h, leaf_temp = Tleaf,
    atm_o2_kpa = 21, atm_kpa = 101.3)
}

test_that("thermal damage is off by default and inert (backward compatible)", {
  l <- mk_leaf(); set_phys(l, 30)
  expect_false(l$use_thermal_damage_)
  expect_equal(l$N_, 1.0)
})

test_that("cold leaf suffers no damage; hot leaf is downscaled by N in (0,1)", {
  cold <- mk_leaf(); cold$use_thermal_damage_ <- TRUE; set_phys(cold, 15)
  expect_equal(cold$N_, 1.0, tolerance = 1e-6)

  hot_off <- mk_leaf(); set_phys(hot_off, 45)
  hot_on  <- mk_leaf(); hot_on$use_thermal_damage_ <- TRUE; set_phys(hot_on, 45)
  expect_true(hot_on$N_ > 0 && hot_on$N_ < 1)
  # jmax is downscaled by exactly N relative to the undamaged leaf
  expect_equal(hot_on$jmax_, hot_off$jmax_ * hot_on$N_, tolerance = 1e-6)
  expect_true(hot_on$electron_transport_ < hot_off$electron_transport_)
})

test_that("N matches the closed-form switch-dominated quasi-steady balance", {
  l <- mk_leaf(); l$use_thermal_damage_ <- TRUE; set_phys(l, 45)
  Tc <- 38; m <- 1; k_d1 <- 864; k_r1 <- 864; m_rep <- 0.4; t_rep <- 45
  Sd <- 1 / (1 + exp(-m * (45 - Tc)))
  Sr <- 1 / (1 + exp(-(-m_rep) * (45 - t_rep)))
  expect_equal(l$N_, (k_r1 * Sr) / (k_r1 * Sr + k_d1 * Sd), tolerance = 1e-9)
})

test_that("repair axis raises N and acclimation input raises T_crit", {
  base <- mk_leaf(); base$use_thermal_damage_ <- TRUE; set_phys(base, 45)
  rep  <- mk_leaf(); rep$use_thermal_damage_ <- TRUE; rep$k_r1_0_ <- 5000; set_phys(rep, 45)
  acc  <- mk_leaf(); acc$use_thermal_damage_ <- TRUE; acc$A_crit_ <- 10; set_phys(acc, 45)
  expect_true(rep$N_ > base$N_)
  expect_true(acc$N_ > base$N_)
})

test_that("tolerance T_opt offset shifts photosynthetic capacity upward", {
  # Strip the damage multiplier to isolate the capacity (d_S) effect.
  l0 <- mk_leaf(); l0$use_thermal_damage_ <- TRUE; set_phys(l0, 40)
  l8 <- mk_leaf(); l8$use_thermal_damage_ <- TRUE; l8$topt_offset_ <- 8; set_phys(l8, 40)
  expect_true((l8$jmax_ / l8$N_) > (l0$jmax_ / l0$N_))
})

test_that("N is finite, bounded and monotone non-increasing across a wide Tleaf sweep", {
  Ts <- seq(-30, 90, by = 1)
  Ns <- vapply(Ts, function(T) { l <- mk_leaf(); l$use_thermal_damage_ <- TRUE; set_phys(l, T); l$N_ }, numeric(1))
  expect_true(all(is.finite(Ns)))
  expect_true(all(Ns >= 0 & Ns <= 1))
  expect_true(all(diff(Ns) <= 1e-9))
})

test_that("outlandish traits and environment still give finite output", {
  l <- mk_leaf(); l$use_thermal_damage_ <- TRUE
  l$topt_offset_ <- 50; l$tcrit_0_ <- 5
  set_phys(l, 99)
  expect_true(is.finite(l$jmax_) && is.finite(l$N_) && is.finite(l$electron_transport_))
})
