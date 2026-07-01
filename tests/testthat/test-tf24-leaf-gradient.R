# TF24 leaf-trait gradient (#472 scope B, Phase F): d(profit*)/d(vcmax_25) at the
# optimised leaf operating point, via forward-mode AD + the implicit function theorem
# at the psi_stem->ci root-find (extending the merged #539 d/d(collar) to a trait) and
# the envelope theorem for the collar optimisation. Compiled into plant.so, exposed on
# the Leaf R class -- runs in plain R, no sourceCpp.

test_that("dprofit_dvcmax25 matches a re-optimising finite difference", {
  rc <- 2.65; rb <- 1.29
  # g1_TF24 = 5 gives a healthy INTERIOR optimum (positive profit, dprofit/dcollar ~ 0),
  # where the envelope theorem applies; the default cost can put the optimum on the
  # feasibility boundary (a constrained optimum) where it does not.
  mk <- function(vc) {
    Leaf(vcmax_25 = vc, jmax_25 = 167, c = 2.04, b = 3, psi_crit = 5,
         root_c = rc, root_b = rb, root_psi_crit = rb * (log(1 / 0.05))^(1 / rc),
         beta2 = 1, hk_s = 75, a = 0.3, curv_fact_elec_trans = 0.7,
         curv_fact_colim = 0.99, GSS_tol_abs = 1e-8, vulnerability_curve_ncontrol = 100,
         ci_abs_tol = 1e-6, ci_niter = 1000, g1_TF24 = 5,
         beta_R_H = 3.4e3, beta_R_V = 9.4e4)
  }
  theta <- 0.000157 * 20; h <- 5
  run <- function(l) {
    l$set_physiology(area_leaf = 0.05, mass_root_prop = 1, rho = 608, a_bio = 0.0245,
                     PPFD = 1500, psi_soil = 0.1, soil_depth = 1,
                     leaf_specific_conductance_max = theta / h, atm_vpd = 1, ca = 40,
                     sapwood_volume_per_leaf_area = theta * h, leaf_temp = 25,
                     atm_o2_kpa = 21, atm_kpa = 101.3)
    l$find_root_collar_psi()
    l
  }
  l0 <- run(mk(100))
  opt <- -l0$root_collar_psi_
  # interior optimum: dprofit/dcollar ~ 0 (envelope theorem applies)
  expect_lt(abs(l0$dprofit_droot_collar_psi(opt)), 1e-3)

  ad <- l0$dprofit_dvcmax25(opt)
  # envelope re-optimising FD: perturb vcmax_25, re-optimise the whole leaf, read profit*.
  profit_at <- function(vc) run(mk(vc))$profit_
  h_fd <- 1e-4 * 100
  fd <- (profit_at(100 + h_fd) - profit_at(100 - h_fd)) / (2 * h_fd)
  expect_equal(ad, fd, tolerance = 1e-4)
})
