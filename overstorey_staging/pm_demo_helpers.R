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
  root_psi_crit <- root_b * (log(1 / 0.05))^(1 / root_c)
  Leaf(vcmax_25 = p$vcmax_25, jmax_25 = p$jmax_25, c = p$c, b = p$b,
       psi_crit = p$psi_crit, root_c = root_c, root_b = root_b,
       root_psi_crit = root_psi_crit, beta2 = p$beta2, a = p$a,
       curv_fact_elec_trans = p$curv_fact_elec_trans,
       curv_fact_colim = p$curv_fact_colim,
       GSS_tol_abs = ctrl$GSS_tol_abs,
       vulnerability_curve_ncontrol = ctrl$vulnerability_curve_ncontrol,
       ci_abs_tol = ctrl$ci_abs_tol, ci_niter = ctrl$ci_niter,
       g1_TF24 = p$g1_TF24, beta_R_H = 3.4e2, beta_R_V = 9.4e3)
}

## Default well-watered, well-rooted, moderate-conductance operating point (so
## the profit optimum opens the stomata; a shut-down leaf makes PM vs Fick
## degenerate). Override via ... to explore.
pm_leaf_config <- function(...) {
  cfg <- list(area_leaf = 1.0, mass_root_prop = 20.0, psi_soil = 0.3,
              leaf_specific_conductance_max = 5e-3, sapwood_volume_per_leaf_area = 1e-3,
              ca = 40, atm_o2_kpa = 21, atm_kpa = 101.3, d = 0.05, wind_speed = 2.0)
  utils::modifyList(cfg, list(...))
}

## Configure physiology on a leaf for one environment cell.
pm_set_physiology <- function(l, PAR, Tair, VPD, pm, cfg = pm_leaf_config()) {
  l$use_energy_balance_ <- isTRUE(pm)
  l$d_ <- cfg$d
  l$wind_speed_ <- cfg$wind_speed
  l$set_physiology(
    area_leaf = cfg$area_leaf, mass_root_prop = cfg$mass_root_prop, rho = 608,
    a_bio = 0.0245, PPFD = PAR, psi_soil = cfg$psi_soil, soil_depth = 1,
    leaf_specific_conductance_max = cfg$leaf_specific_conductance_max,
    atm_vpd = VPD, ca = cfg$ca,
    sapwood_volume_per_leaf_area = cfg$sapwood_volume_per_leaf_area,
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
             opt_psi_stem = l$opt_psi_stem_, A = l$assim_colimited_,
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

## Profit "anatomy": scan profit / assimilation / hydraulic cost across candidate
## psi_stem at one environment, for Fick and PM. psi_upstream is fixed at the
## configured soil potential (single-layer supply).
pm_profit_curve <- function(PAR, Tair, VPD, pm, psi_stem_seq, cfg = pm_leaf_config()) {
  l <- pm_make_leaf()
  pm_set_physiology(l, PAR, Tair, VPD, pm, cfg)
  psi_up <- cfg$psi_soil
  rows <- lapply(psi_stem_seq, function(ps) {
    profit <- l$profit_psi_stem_TF(ps, psi_up)  # sets leaf states + cost as a side effect
    data.frame(pm = pm, psi_stem = ps, profit = profit,
               assim = l$assim_colimited_, hydraulic_cost = l$hydraulic_cost_)
  })
  do.call(rbind, rows)
}
