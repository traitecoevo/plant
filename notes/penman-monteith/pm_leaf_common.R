## Shared leaf driver for the Penman-Monteith evaluation scripts in this folder
## (step5_fick_vs_pm.R, sensitivity_leaf.R). One definition of the leaf, so the
## Step 5 re-run and the 7.5.2/7.5.3 sensitivity tests solve exactly the same
## leaf and differ only in the perturbation they inject.
##
## Run every script from the plant repo root with the pinned stack installed
## (see pm_check_stack); plant is loaded with library(), not load_all().

suppressPackageStartupMessages(library(plant))

## Print and check the stack: plant's location and the three versions, against
## the `==` pins in the repo's DESCRIPTION LinkingTo when it is readable.
pm_check_stack <- function(description = "DESCRIPTION") {
  cat("plant loaded from:", find.package("plant"), "\n")
  v <- vapply(c("plant", "phylloptim", "odelia"),
              function(p) as.character(utils::packageVersion(p)), "")
  cat(sprintf("%s %s\n", names(v), v), sep = "")
  if (file.exists(description)) {
    d <- read.dcf(description, fields = c("Version", "LinkingTo"))
    stopifnot(v[["plant"]] == unname(d[1, "Version"]))
    for (p in c("phylloptim", "odelia")) {
      pin <- regmatches(d[1, "LinkingTo"],
                        regexpr(paste0(p, " \\(== [0-9.]+\\)"), d[1, "LinkingTo"]))
      want <- unname(sub(".*== ([0-9.]+)\\)", "\\1", pin))
      if (length(want) == 1L && v[[p]] != want) {
        stop(p, " ", v[[p]], " installed but DESCRIPTION pins ", want)
      }
    }
  }
  invisible(v)
}

## PM constants, mirroring phylloptim/constants.hpp. Used only to read the
## operating-point leaf temperature back out in R: plant's Leaf binding exposes
## Tair_, Rn_, ra_ and transpiration_ but not Tleaf_ itself.
PM_LAMBDA   <- 2.45e6   # latent heat of vaporisation, J kg^-1
PM_RHO_CP   <- 1200.0   # volumetric heat capacity of air, J m^-3 K^-1
PM_T_MIN    <- -40      # leaf_temp_min, deg C
PM_T_MAX    <- 70       # leaf_temp_max, deg C
PM_VPD_MIN  <- 0.01     # vpd_leaf_min, kPa
PM_C_RA     <- 200      # aerodynamic_resistance_coef, s^0.5 m^-1
PM_RA_FIXED <- 50       # aerodynamic_resistance_fixed, s m^-1

pm_esat <- function(T) 0.6108 * exp(17.27 * T / (T + 237.3))   # kPa, Tetens

## The TF24-default strategy the leaf is built from.
PM_PARS <- TF24_Strategy()$pars
PM_CTRL <- Control()

## A Leaf exactly as TF24_Strategy::prepare_strategy builds it: same
## constructor slots, R_d_25 at its default (asserted in pm_solve_cell), and
## the TF24_floor curve at lambda_o = 0 on the collar route.
pm_make_leaf <- function(p = PM_PARS, ctrl = PM_CTRL) {
  l <- Leaf(vcmax_25 = p$vcmax_25, stem_c = p$stem_c, stem_P50 = p$stem_P50,
            root_c = p$root_c, root_P50 = p$root_P50,
            TF24_beta2 = p$TF24_beta2, jmax_25 = p$jmax_25, a = p$a,
            curv_fact_elec_trans = p$curv_fact_elec_trans,
            curv_fact_colim = p$curv_fact_colim,
            GSS_tol_abs = ctrl$GSS_tol_abs,
            vulnerability_curve_ncontrol = ctrl$vulnerability_curve_ncontrol,
            ci_abs_tol = ctrl$ci_abs_tol, ci_niter = ctrl$ci_niter,
            TF24_cost_scale = p$TF24_cost_scale)
  l$TF24_floor_lambda_o <- p$TF24_floor_lambda_o
  l$set_model("TF24_floor", "collar")
  l
}

## The Step 5 single-leaf operating point: well watered, well rooted, moderate
## conductance, so the optimum opens the stomata. July's script also passed
## sapwood_volume_per_leaf_area = 1e-3; set_physiology no longer takes it.
PM_GEOM <- list(root_carbon_per_leaf_area = 20, psi_soil = 0.3,
                leaf_specific_conductance_max = 5e-3, ca = 40,
                atm_o2_kpa = 21, atm_kpa = 101.3)

## A leaf with physiology set for one cell. `ra` and `rn_scale` are the
## sensitivity injections: they overwrite the boundary-layer resistance and
## scale the net radiation AFTER set_physiology has derived them and BEFORE the
## solve reads them (both are read only inside the solve, by leaf_temp_from_E).
pm_cell_leaf <- function(PAR, Tair, VPD, pm, d = 0.05, U = 2.0,
                         ra = NULL, rn_scale = 1, geom = PM_GEOM) {
  l <- pm_make_leaf()
  l$use_energy_balance_ <- pm
  l$d_ <- d
  l$wind_speed_ <- U
  stopifnot(identical(l$use_energy_balance_, pm))
  l$set_physiology(
    root_network = phylloptim::root_network_from_carbon(
      geom$root_carbon_per_leaf_area, soil_depth = 1,
      beta_R_H = 3.4e2, beta_R_V = 9.4e3),
    PPFD = PAR, psi_soil = geom$psi_soil, soil_depth = 1,
    leaf_specific_conductance_max = geom$leaf_specific_conductance_max,
    atm_vpd = VPD, ca = geom$ca, leaf_temp = Tair,
    atm_o2_kpa = geom$atm_o2_kpa, atm_kpa = geom$atm_kpa)
  attr(l, "ra_model") <- l$ra_
  attr(l, "Rn_model") <- l$Rn_
  if (!is.null(ra)) l$ra_ <- ra
  if (rn_scale != 1) l$Rn_ <- rn_scale * l$Rn_
  l
}

## The profit maximum over the whole collar bracket, by brute force: a scan of
## evaluate_root_collar_psi() over [psi_soil, root_psi_crit] and a golden-
## section refinement around the best scan point. Independent of the solver, so
## it says whether the solver's answer is the GLOBAL maximum. On the PM path the
## profit can have two maxima -- transpire hard to cool, or shut -- and the
## solver does not always return the higher one (see the 2026-10 write-up).
pm_scan_collar <- function(l, geom = PM_GEOM, n = 300) {
  f <- function(x) tryCatch(l$evaluate_root_collar_psi(x),
                            error = function(e) -Inf)
  xs <- seq(geom$psi_soil + 1e-5, PM_PARS$root_psi_crit, length.out = n)
  ps <- vapply(xs, f, 0)
  j <- which.max(ps)
  o <- stats::optimize(f, c(xs[max(j - 1L, 1L)], xs[min(j + 1L, n)]),
                       maximum = TRUE, tol = 1e-9)
  if (o$objective >= ps[j]) c(collar = o$maximum, profit = o$objective)
  else c(collar = xs[j], profit = ps[j])
}

## Solve one cell with plant's solver, as TF24 does. With `global = TRUE`, also
## run pm_scan_collar and, where it beats the solver's profit by more than
## `tol`, report the leaf at the scan's maximum instead: `gap` is that profit
## difference (umol m-2 s-1) and `solver_profit` the solver's own value.
pm_solve_cell <- function(PAR, Tair, VPD, pm, d = 0.05, U = 2.0,
                          ra = NULL, rn_scale = 1, geom = PM_GEOM,
                          global = FALSE, tol = 1e-4) {
  mk <- function() pm_cell_leaf(PAR, Tair, VPD, pm, d, U, ra, rn_scale, geom)
  l <- mk()
  ra_model <- attr(l, "ra_model"); Rn_model <- attr(l, "Rn_model")
  l$optimise()
  solver_profit <- l$profit_
  gap <- NA_real_
  if (global) {
    l2 <- mk()
    s <- pm_scan_collar(l2, geom)
    gap <- max(s[["profit"]] - solver_profit, 0)
    if (gap > tol) {
      l2$evaluate_root_collar_psi(s[["collar"]])
      l <- l2
    }
  }
  E <- l$transpiration_
  Tleaf <- if (pm) {
    min(max(l$Tair_ + (l$Rn_ - PM_LAMBDA * E) * l$ra_ / PM_RHO_CP, PM_T_MIN), PM_T_MAX)
  } else Tair
  vpd_leaf <- if (pm) max(VPD + pm_esat(Tleaf) - pm_esat(Tair), PM_VPD_MIN) else VPD
  data.frame(PAR = PAR, Tair = Tair, VPD = VPD, pm = pm,
             Rn = l$Rn_, ra = l$ra_, ra_model = ra_model, Rn_model = Rn_model,
             Tleaf = Tleaf, vpd_leaf = vpd_leaf,
             E = E, A = l$assim_colimited_, gs = l$stom_cond_CO2_,
             R_d = l$R_d_, ci = l$ci_, opt_psi_stem = l$opt_psi_stem_,
             profit = if (global && gap > tol) solver_profit + gap else l$profit_,
             solver_profit = solver_profit, gap = gap)
}
