# Leaf thermal-damage / acclimation layer (ATLS merged revision; issue #566).
# The merged revision unifies the reversible high-T decline onto the Medlyn
# peaked-Arrhenius curve (the deactivated fraction phi_d = K/(1+K) comes from the
# SAME K(T) already in that curve) and drives lasting damage through two pools
# I_r/I_p that the strategy integrates. At the Leaf level we test: the emergent
# damage onset, the phi_d curve, the (1-I_r-I_p) discount on BOTH capacities, the
# repair-collapse gate k_rec_eff, the numerical-stability contract, and the
# merged R_d_ cost terms.

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

# Closed-form deactivated fraction for the default (unshifted) damage curve:
# K = exp((d_S*T_K - H_d)/(R*T_K)), phi_d = K/(1+K), with d_S = 650, H_d = 200000.
phi_d_ref <- function(Tleaf, d_S = 650, H_d = 200000, R = 8.314, T0 = 273.15) {
  T_K <- Tleaf + T0
  K <- exp((d_S * T_K - H_d) / (R * T_K))
  K / (1 + K)
}

test_that("thermal damage is off by default and inert (backward compatible)", {
  l <- mk_leaf(); set_phys(l, 30)
  expect_false(l$use_thermal_damage_)
  expect_equal(l$phi_d_, 0.0)     # no damage machinery evaluated
  # Gate on but undamaged (I_r=I_p=0): capacity is identical to the gate-off leaf
  # (lasting damage is a state, so a fresh hot leaf loses nothing instantly beyond
  # the reversible peaked-Arrhenius response already in f_peak).
  off <- mk_leaf(); set_phys(off, 45)
  on  <- mk_leaf(); on$use_thermal_damage_ <- TRUE; set_phys(on, 45)
  expect_equal(on$vcmax_, off$vcmax_, tolerance = 1e-9)
  expect_equal(on$jmax_,  off$jmax_,  tolerance = 1e-9)
})

test_that("phi_d matches the Medlyn deactivation equilibrium; onset is emergent", {
  # phi_d = K/(1+K) from the same K(T) in peak_arrh_curve; no tcrit_0 parameter.
  for (T in c(20, 34.54, 45, 55)) {
    l <- mk_leaf(); l$use_thermal_damage_ <- TRUE; set_phys(l, T)
    expect_equal(l$phi_d_, phi_d_ref(T), tolerance = 1e-9)
  }
  # Emergent half-unfolding (K=1) sits at T_K = H_d/d_S = 200000/650 ~ 34.54 C.
  l <- mk_leaf(); l$use_thermal_damage_ <- TRUE; set_phys(l, 200000 / 650 - 273.15)
  expect_equal(l$phi_d_, 0.5, tolerance = 1e-6)
})

test_that("lasting damage discounts BOTH Vcmax and Jmax by (1 - I_r - I_p)", {
  undmg <- mk_leaf(); undmg$use_thermal_damage_ <- TRUE; set_phys(undmg, 30)
  dmg   <- mk_leaf(); dmg$use_thermal_damage_ <- TRUE
  dmg$I_r_ <- 0.3; dmg$I_p_ <- 0.1                 # surv = 0.6
  set_phys(dmg, 30)
  expect_equal(dmg$vcmax_ / undmg$vcmax_, 0.6, tolerance = 1e-9)
  expect_equal(dmg$jmax_  / undmg$jmax_,  0.6, tolerance = 1e-9)
  # Bounded: I_r + I_p > 1 floors the survival fraction at a tiny positive value
  # (strictly positive capacity keeps the photosynthesis solve well-conditioned).
  crush <- mk_leaf(); crush$use_thermal_damage_ <- TRUE
  crush$I_r_ <- 0.8; crush$I_p_ <- 0.5; set_phys(crush, 30)
  expect_equal(crush$vcmax_ / undmg$vcmax_, 1e-6, tolerance = 1e-9)
  expect_equal(crush$jmax_  / undmg$jmax_,  1e-6, tolerance = 1e-9)
})

test_that("thermostability offset shifts the whole response up in temperature", {
  # topt_offset raises T_opt (capacity) AND the damage onset together: at a fixed
  # hot leaf, a more thermostable leaf has a lower deactivated fraction.
  l0 <- mk_leaf(); l0$use_thermal_damage_ <- TRUE; set_phys(l0, 40)
  l8 <- mk_leaf(); l8$use_thermal_damage_ <- TRUE; l8$topt_offset_ <- 8; set_phys(l8, 40)
  expect_gt(l8$vcmax_, l0$vcmax_)   # capacity optimum moved up
  expect_lt(l8$phi_d_, l0$phi_d_)   # damage onset moved up (less unfolded at 40 C)
})

test_that("restorative rate k_rec_eff is gated off above t_rep_cut", {
  cold <- mk_leaf(); cold$use_thermal_damage_ <- TRUE; set_phys(cold, 20)
  hot  <- mk_leaf(); hot$use_thermal_damage_ <- TRUE; set_phys(hot, 60)
  expect_equal(cold$k_rec_eff_, cold$k_rec_, tolerance = 1e-3)  # background rate when cool
  expect_lt(hot$k_rec_eff_, 1e-2 * hot$k_rec_)                  # collapses when too hot
  # At exactly t_rep_cut the logistic is 0.5, so k_rec_eff = k_rec/2.
  cut <- mk_leaf(); cut$use_thermal_damage_ <- TRUE; set_phys(cut, cut$t_rep_cut_)
  expect_equal(cut$k_rec_eff_, cut$k_rec_ / 2, tolerance = 1e-9)
})

test_that("phi_d is finite, bounded [0,1] and monotone increasing across a wide sweep", {
  Ts <- seq(-30, 90, by = 1)
  phis <- vapply(Ts, function(T) {
    l <- mk_leaf(); l$use_thermal_damage_ <- TRUE; set_phys(l, T); l$phi_d_
  }, numeric(1))
  expect_true(all(is.finite(phis)))
  expect_true(all(phis >= 0 & phis <= 1))
  expect_true(all(diff(phis) >= -1e-9))   # monotone non-decreasing in temperature
})

test_that("outlandish traits and environment still give finite output", {
  l <- mk_leaf(); l$use_thermal_damage_ <- TRUE
  l$topt_offset_ <- 50; l$I_r_ <- 0.4; l$I_p_ <- 0.4
  set_phys(l, 99)
  expect_true(is.finite(l$vcmax_) && is.finite(l$jmax_) &&
              is.finite(l$phi_d_) && is.finite(l$electron_transport_))
  expect_true(l$phi_d_ >= 0 && l$phi_d_ <= 1)
})

# Midday evaluation. On the Penman-Monteith path the damage substrate phi_d is
# recomputed at the operating-point leaf temperature Tleaf = f(E) inside
# set_leaf_states_rates_from_psi_stem, so transpirational cooling lowers phi_d
# (and hence the I_r leak the strategy integrates). This is the "avoidance" axis:
# investing in transpiration lowers midday Tleaf and reduces the damage flux,
# self-selecting inside the psi_stem optimiser.
test_that("on the PM path phi_d is evaluated at the operating-point Tleaf (avoidance)", {
  mk_pm <- function() {
    l <- mk_leaf()
    l$use_energy_balance_ <- TRUE
    l$use_thermal_damage_ <- TRUE
    l$d_ <- 0.05; l$wind_speed_ <- 2.0
    set_phys(l, 42)   # hot, bright midday: damage substrate active and E-sensitive
    l
  }
  # More transpiration => cooler operating-point leaf => lower phi_d.
  lo <- mk_pm(); lo$set_leaf_states_rates_from_psi_stem(1.0, 0.2)
  hi <- mk_pm(); hi$set_leaf_states_rates_from_psi_stem(3.5, 0.2)
  expect_gt(hi$transpiration_, lo$transpiration_)   # bigger psi_stem draw -> more E
  expect_lt(hi$phi_d_, lo$phi_d_)                    # ... which cools & lowers phi_d
  expect_true(hi$phi_d_ >= 0 && hi$phi_d_ <= 1 && is.finite(hi$phi_d_))
  # Sanity: with damage off, phi_d stays 0 regardless of the operating point.
  off <- mk_leaf(); off$use_energy_balance_ <- TRUE; set_phys(off, 42)
  off$set_leaf_states_rates_from_psi_stem(3.5, 0.2)
  expect_equal(off$phi_d_, 0.0)
})

# The merged thermal maintenance/activity respiration costs added to R_d_.
# Coefficients default 0, so a bare Leaf pays nothing (R_d_ = vcmax_*0.015);
# TF24t sets them. Each term maps to one axis of the merged cost structure.
test_that("thermal maintenance/activity costs raise R_d_ (and default to inert)", {
  # Default coefficients (0) -> R_d_ is exactly the base dark respiration.
  base <- mk_leaf(); base$use_thermal_damage_ <- TRUE; set_phys(base, 30)
  expect_equal(base$R_d_, base$vcmax_ * 0.015, tolerance = 1e-9)

  # Acclimation maintenance ~ held thermostability load A.
  acc <- mk_leaf(); acc$use_thermal_damage_ <- TRUE
  acc$A_ <- 3; acc$c_acclim_maint_ <- 0.1
  set_phys(acc, 30)
  expect_equal(acc$R_d_ - acc$vcmax_ * 0.015, 0.1 * 3, tolerance = 1e-9)

  # Repair standing maintenance ~ resynthesis capacity k_rec (paid even when cold).
  rm <- mk_leaf(); rm$use_thermal_damage_ <- TRUE
  rm$c_repair_maint_ <- 1e-2
  set_phys(rm, 10)
  expect_equal(rm$R_d_ - rm$vcmax_ * 0.015, 1e-2 * rm$k_rec_, tolerance = 1e-9)

  # Protection standing maintenance ~ (k_i_ref - k_i), paid for buying protection.
  pr <- mk_leaf(); pr$use_thermal_damage_ <- TRUE
  pr$k_i_ <- 0.2; pr$c_protect_maint_ <- 0.5   # k_i_ref default 0.5 -> protection 0.3
  set_phys(pr, 10)
  expect_equal(pr$R_d_ - pr$vcmax_ * 0.015, 0.5 * (pr$k_i_ref_ - 0.2), tolerance = 1e-9)

  # Repair activity ~ realized resynthesis flux k_rec_eff * I_r, only when damaged.
  ra <- mk_leaf(); ra$use_thermal_damage_ <- TRUE
  ra$c_repair_flux_ <- 1.0; ra$I_r_ <- 0.4
  set_phys(ra, 30)
  expect_equal(ra$R_d_ - ra$vcmax_ * 0.015, 1.0 * ra$k_rec_eff_ * 0.4, tolerance = 1e-9)
  expect_gt(ra$R_d_, ra$vcmax_ * 0.015)  # activity cost genuinely nonzero when damaged
})
