## Step 5 (doc 7.5.1): factorial Fick vs Penman-Monteith at the leaf level.
## Drives a single TF24-default leaf across Tair x VPD x PAR, PM off vs on, and
## reports dTleaf, dE, dA, dgs. Writes a scorecard RDS + CSV.
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
outstem <- if (length(args) >= 1) args[[1]] else "step5"

## PM physical constants (mirror leaf_model.h) for the R-side Tleaf readout.
LAMBDA <- 2.45e6      # J kg^-1
RHO_CP <- 1200.0      # J m^-3 K^-1

## Build a Leaf with TF24 defaults (so it actually transpires), matching how
## TF24_Strategy::prepare_strategy constructs it.
s <- TF24_Strategy(); p <- s$pars
ctrl <- Control()
root_c <- 2.680147; root_b <- 3.898245
root_psi_crit <- root_b * (log(1 / 0.05))^(1 / root_c)
mk_leaf <- function() {
  Leaf(vcmax_25 = p$vcmax_25, jmax_25 = p$jmax_25, c = p$c, b = p$b,
       psi_crit = p$psi_crit, root_c = root_c, root_b = root_b,
       root_psi_crit = root_psi_crit, beta2 = p$beta2, a = p$a,
       curv_fact_elec_trans = p$curv_fact_elec_trans,
       curv_fact_colim = p$curv_fact_colim,
       GSS_tol_abs = ctrl$GSS_tol_abs,
       vulnerability_curve_ncontrol = ctrl$vulnerability_curve_ncontrol,
       ci_abs_tol = ctrl$ci_abs_tol, ci_niter = ctrl$ci_niter,
       g1_TF24 = p$g1_TF24)
}

## Representative single-leaf geometry, tuned so the profit optimum opens the
## stomata (a shut-down leaf makes the Fick-vs-PM contrast degenerate). These are
## a well-watered, well-rooted, moderate-conductance operating point.
lsc_max <- 5e-3; svla <- 1e-3
area_leaf <- 1.0; mass_root <- 20.0; psi_soil <- 0.3

solve_cell <- function(PAR, Tair, VPD, pm, d = 0.05, U = 2.0) {
  l <- mk_leaf()
  l$use_energy_balance_ <- pm; l$d_ <- d; l$wind_speed_ <- U
  l$set_physiology(area_leaf = area_leaf, mass_root_prop = mass_root, rho = p$rho,
                   a_bio = p$a_bio, PPFD = PAR, psi_soil = psi_soil, soil_depth = 1,
                   leaf_specific_conductance_max = lsc_max, atm_vpd = VPD, ca = 40,
                   sapwood_volume_per_leaf_area = svla, leaf_temp = Tair,
                   atm_o2_kpa = 21, atm_kpa = 101.3)
  l$find_root_collar_psi()
  Tleaf <- if (pm) l$Tair_ + (l$Rn_ - LAMBDA * l$transpiration_) * l$ra_ / RHO_CP else Tair
  list(E = l$transpiration_, A = l$assim_colimited_, gs = l$stom_cond_CO2_,
       Tleaf = Tleaf, opt_psi = l$opt_psi_stem_, profit = l$profit_)
}

grid <- expand.grid(PAR = c(200, 500, 1000, 2000),
                    Tair = c(15, 25, 35, 40),
                    VPD = c(0.5, 1, 2, 3))
rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  off <- solve_cell(g$PAR, g$Tair, g$VPD, FALSE)
  on  <- solve_cell(g$PAR, g$Tair, g$VPD, TRUE)
  # Relative dA only where the Fick benchmark A is materially positive (>1 umol);
  # near the compensation point the ratio is a meaningless artifact, so report
  # absolute dA there instead.
  relA <- if (is.finite(off$A) && off$A > 1.0) 100 * (on$A - off$A) / off$A else NA_real_
  rows[[i]] <- data.frame(
    PAR = g$PAR, Tair = g$Tair, VPD = g$VPD,
    Tleaf_pm = on$Tleaf, dTleaf = on$Tleaf - g$Tair,
    E_off = off$E, E_on = on$E, dE = on$E - off$E,
    A_off = off$A, A_on = on$A, dA = on$A - off$A, relA_pct = relA,
    gs_off = off$gs, gs_on = on$gs, dgs = on$gs - off$gs,
    profit_off = off$profit, profit_on = on$profit)
}
res <- do.call(rbind, rows)

saveRDS(res, paste0(outstem, ".rds"))
write.csv(res, paste0(outstem, ".csv"), row.names = FALSE)

## Gate summary (doc 7.4 / 7.5): keep PM if max|dTleaf| > 2C AND max|relA| > 5%.
fin <- res[is.finite(res$relA_pct), ]
cat(sprintf("cells: %d (%d with finite relA)\n", nrow(res), nrow(fin)))
cat(sprintf("max |Tleaf-Tair|      = %.2f C\n", max(abs(res$dTleaf), na.rm = TRUE)))
cat(sprintf("max |relative dA|     = %.1f %%\n", max(abs(fin$relA_pct), na.rm = TRUE)))
cat(sprintf("median |Tleaf-Tair|   = %.2f C\n", median(abs(res$dTleaf), na.rm = TRUE)))
cat("--- hot/high-PAR/high-VPD corner (Tair>=35, PAR>=1000, VPD>=2) ---\n")
corner <- subset(res, Tair >= 35 & PAR >= 1000 & VPD >= 2)
print(round(corner[, c("PAR","Tair","VPD","dTleaf","relA_pct","dE","dgs")], 3), row.names = FALSE)
cat("--- cool/low-PAR corner (Tair<=25, PAR<=500) mean |dTleaf| ---\n")
cool <- subset(res, Tair <= 25 & PAR <= 500)
cat(sprintf("mean |dTleaf| = %.3f C, mean |relA| = %.2f %%\n",
            mean(abs(cool$dTleaf)), mean(abs(cool$relA_pct[is.finite(cool$relA_pct)]))))
