## Shared leaf-driving helpers for the TF24 Penman-Monteith leaf demo (#523).
##
## Kept separate from the .qmd so the smoke test (tests/testthat/test-pm-leaf-demo.R)
## can exercise exactly the same leaf-driving code as the vignette, guaranteeing
## the demo cannot silently rot. No plotting here — pure model driving.
##
## Assumes the `plant` package is loaded (library(plant) or devtools::load_all()).

## PM leaf-energy-balance physical constants (mirror inst/include/plant/leaf_model.h),
## used only to recover the operating-point leaf temperature Tleaf in R.
pm_lambda <- 2.45e6   # latent heat of vaporisation, J kg^-1
pm_rho_cp <- 1200.0   # volumetric heat capacity of air, J m^-3 K^-1

## Build a Leaf with TF24 defaults, matching TF24_Strategy::prepare_strategy.
pm_make_leaf <- function() {
  s <- TF24_Strategy(); p <- s$pars; ctrl <- Control()
  root_c <- 2.680147; root_b <- 3.898245
  # phylloptim 0.6.0 parameterises both curves on P50 and derives stem_b and psi_crit
  # itself; the root pair is a literal Weibull scale (root_b) here, so convert
  # (#622). Argument names match phylloptim's since #634.
  root_P50 <- root_b * (log(2))^(1 / root_c)
  Leaf(vcmax_25 = p$vcmax_25, jmax_25 = p$jmax_25, stem_c = p$stem_c, stem_P50 = p$stem_P50,
       root_c = root_c, root_P50 = root_P50, TF24_beta2 = p$TF24_beta2, a = p$a,
       curv_fact_elec_trans = p$curv_fact_elec_trans,
       curv_fact_colim = p$curv_fact_colim,
       GSS_tol_abs = ctrl$GSS_tol_abs,
       vulnerability_curve_ncontrol = ctrl$vulnerability_curve_ncontrol,
       ci_abs_tol = ctrl$ci_abs_tol, ci_niter = ctrl$ci_niter,
       TF24_cost_scale = p$TF24_cost_scale)
}

## Default well-watered, well-rooted, moderate-conductance operating point (so
## the profit optimum opens the stomata; a shut-down leaf makes PM vs Fick
## degenerate). Override via ... to explore.
pm_leaf_config <- function(...) {
  cfg <- list(root_carbon_per_leaf_area = 20.0, psi_soil = 0.3,
              leaf_specific_conductance_max = 5e-3,
              ca = 40, atm_o2_kpa = 21, atm_kpa = 101.3, d = 0.05, wind_speed = 2.0)
  utils::modifyList(cfg, list(...))
}

## Configure physiology on a leaf for one environment cell.
pm_set_physiology <- function(l, PAR, Tair, VPD, pm, cfg = pm_leaf_config()) {
  l$use_energy_balance_ <- isTRUE(pm)
  l$d_ <- cfg$d
  l$wind_speed_ <- cfg$wind_speed
  # phylloptim #33: the leaf takes the resistances, so the architecture model runs
  # here. Same constants TF24_Strategy uses.
  l$set_physiology(
    root_network = phylloptim::root_network_from_carbon(
      cfg$root_carbon_per_leaf_area, soil_depth = 1,
      beta_R_H = 3.4e2, beta_R_V = 9.4e3),
    PPFD = PAR, psi_soil = cfg$psi_soil, soil_depth = 1,
    leaf_specific_conductance_max = cfg$leaf_specific_conductance_max,
    atm_vpd = VPD, ca = cfg$ca,
    leaf_temp = Tair, atm_o2_kpa = cfg$atm_o2_kpa, atm_kpa = cfg$atm_kpa)
  invisible(l)
}

## Operating-point leaf temperature (deg C). On the PM path, recovered from the
## exposed energy-balance fields; on the Fick path Tleaf == Tair by definition.
pm_leaf_temp <- function(l, Tair) {
  if (isTRUE(l$use_energy_balance_)) {
    l$Tair_ + (l$Rn_ - pm_lambda * l$transpiration_) * l$ra_ / pm_rho_cp
  } else {
    Tair
  }
}

## Solve the optimal operating point at one environment cell; returns a 1-row
## data.frame of leaf outputs.
pm_solve_cell <- function(PAR, Tair, VPD, pm, cfg = pm_leaf_config()) {
  l <- pm_make_leaf()
  pm_set_physiology(l, PAR, Tair, VPD, pm, cfg)
  l$find_root_collar_psi()
  data.frame(PAR = PAR, Tair = Tair, VPD = VPD, pm = pm,
             Tleaf = pm_leaf_temp(l, Tair),
             opt_psi_stem = l$opt_psi_stem_,
             opt_root_psi = l$opt_root_psi_,
             A = l$assim_colimited_,
             gs = l$stom_cond_CO2_, E = l$transpiration_, profit = l$profit_)
}

## Solve a whole environment grid, Fick and PM, returning a long data.frame.
pm_solve_grid <- function(grid, cfg = pm_leaf_config()) {
  out <- vector("list", nrow(grid) * 2L); k <- 0L
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    for (pm in c(FALSE, TRUE)) {
      k <- k + 1L
      out[[k]] <- pm_solve_cell(g$PAR, g$Tair, g$VPD, pm, cfg)
    }
  }
  do.call(rbind, out)
}

## Profit "anatomy" driven by root-collar water potential -- the actual GSS
## search variable inside find_root_collar_psi(). evaluate_root_collar_psi()
## derives psi_stem from each candidate collar potential via the same
## hydraulic transport (find_psi_stem_from_psi_root) the solver itself uses,
## then evaluates profit there -- this is how the model is actually driven.
## (An earlier version of this helper, pm_profit_curve, scanned psi_stem
## directly with psi_upstream fixed at psi_soil; that only agrees with the
## solver when root resistance happens to be negligible, and was removed.)
##
## root_psi_seq is a POSITIVE magnitude -- as is every psi in the leaf package now
## (phylloptim #25), including the reported opt_root_psi; values outside the feasible
## interval for this leaf/environment are
## clamped by evaluate_root_collar_psi rather than extrapolated.
pm_collar_curve <- function(PAR, Tair, VPD, pm, root_psi_seq, cfg = pm_leaf_config()) {
  l <- pm_make_leaf()
  pm_set_physiology(l, PAR, Tair, VPD, pm, cfg)
  rows <- lapply(root_psi_seq, function(rp) {
    profit <- l$evaluate_root_collar_psi(rp)
    data.frame(pm = pm, root_psi = rp, psi_stem = l$opt_psi_stem_,
               profit = profit, assim = l$assim_colimited_,
               hydraulic_cost = l$hydraulic_cost_, E = l$transpiration_,
               leaf_temp = pm_leaf_temp(l, Tair))
  })
  do.call(rbind, rows)
}
