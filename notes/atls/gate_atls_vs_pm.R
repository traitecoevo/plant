## ATLS verification gate (#566, mirrors the #523 Step-5 factorial): does the
## Lumry-Eyring leaf thermal-damage feedback materially change carbon gain vs
## PM-only, under representative Australian heatwave leaf temperatures? Four
## self-contained measurements:
##   (A) CLEAN mechanism isolation -- prescribed leaf temperature (Fick path,
##       Tleaf = driver, no PM overheating confound), damage OFF vs ON, sweeping
##       Tleaf across mild -> heatwave. This is the pure damage response.
##   (B) REALISTIC PM operating point -- use_energy_balance on, N evaluated at the
##       midday operating-point Tleaf (Phase 4 wiring), over Tair x PAR x VPD.
##       Honest about the minimal-cut PM's leaf-overheating behaviour.
##   (C) WHOLE-PLANT net production across the Tair gradient: TF24 (PM-only) vs
##       TF24t-damage-only (costs zeroed, isolating jmax*N) vs full TF24t.
##   (D) ACCLIMATION buffering: does raising the T_crit acclimation state recover
##       N/A at a heatwave leaf temperature within acclimation's reach?
## Writes gate_atls_vs_pm.rds + .csv and prints a keep/drop scorecard.
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
outstem <- if (length(args) >= 1) args[[1]] else "gate_atls_vs_pm"

LAMBDA <- 2.45e6      # J kg^-1   (mirror leaf_model.h, for the R-side Tleaf readout)
RHO_CP <- 1200.0      # J m^-3 K^-1

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
       g1_TF24 = p$g1_TF24, beta_R_H = 3.4e2, beta_R_V = 9.4e3)
}
lsc_max <- 5e-3; svla <- 1e-3
area_leaf <- 1.0; mass_root <- 20.0; psi_soil <- 0.3

## One leaf solve. pm = FALSE prescribes Tleaf = leaf_temp (clean isolation);
## pm = TRUE runs the energy balance and evaluates N at the operating point.
solve_cell <- function(PAR, leaf_temp, VPD, atls, pm, A_crit = 0,
                       d = 0.05, U = 2.0) {
  l <- mk_leaf()
  l$use_energy_balance_ <- pm
  l$use_thermal_damage_ <- atls
  l$A_crit_ <- A_crit
  l$d_ <- d; l$wind_speed_ <- U
  l$set_physiology(area_leaf = area_leaf, mass_root_prop = mass_root, rho = p$rho,
                   a_bio = p$a_bio, PPFD = PAR, psi_soil = psi_soil, soil_depth = 1,
                   leaf_specific_conductance_max = lsc_max, atm_vpd = VPD, ca = 40,
                   sapwood_volume_per_leaf_area = svla, leaf_temp = leaf_temp,
                   atm_o2_kpa = 21, atm_kpa = 101.3)
  l$find_root_collar_psi()
  Tleaf <- if (pm) l$Tair_ + (l$Rn_ - LAMBDA * l$transpiration_) * l$ra_ / RHO_CP
           else leaf_temp
  list(E = l$transpiration_, A = l$assim_colimited_, N = l$N_,
       Tleaf = Tleaf, profit = l$profit_)
}

## ---------------------------------------------------------------------------
## (A) Clean mechanism isolation: prescribed Tleaf, damage OFF vs ON.
gridA <- expand.grid(PAR = c(1000, 2000),
                     Tleaf = c(25, 30, 34, 36, 38, 40, 42, 44, 46),
                     VPD = 2)
rowsA <- lapply(seq_len(nrow(gridA)), function(i) {
  g <- gridA[i, ]
  off <- solve_cell(g$PAR, g$Tleaf, g$VPD, atls = FALSE, pm = FALSE)
  on  <- solve_cell(g$PAR, g$Tleaf, g$VPD, atls = TRUE,  pm = FALSE)
  relA <- if (is.finite(off$A) && off$A > 1) 100 * (on$A - off$A) / off$A else NA_real_
  data.frame(PAR = g$PAR, Tleaf = g$Tleaf, N = on$N,
             A_off = off$A, A_on = on$A, relA_pct = relA)
})
resA <- do.call(rbind, rowsA)

## ---------------------------------------------------------------------------
## (B) Realistic PM operating point: N at the midday operating-point Tleaf.
gridB <- expand.grid(PAR = c(1000, 2000), Tair = c(25, 32, 38, 42, 46), VPD = 2)
rowsB <- lapply(seq_len(nrow(gridB)), function(i) {
  g <- gridB[i, ]
  off <- solve_cell(g$PAR, g$Tair, g$VPD, atls = FALSE, pm = TRUE)
  on  <- solve_cell(g$PAR, g$Tair, g$VPD, atls = TRUE,  pm = TRUE)
  relA <- if (is.finite(off$A) && off$A > 1) 100 * (on$A - off$A) / off$A else NA_real_
  data.frame(PAR = g$PAR, Tair = g$Tair, Tleaf_op = on$Tleaf, N = on$N,
             A_off = off$A, A_on = on$A, relA_pct = relA)
})
resB <- do.call(rbind, rowsB)

saveRDS(list(A = resA, B = resB), paste0(outstem, ".rds"))
write.csv(resA, paste0(outstem, "_A.csv"), row.names = FALSE)
write.csv(resB, paste0(outstem, "_B.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## (C) Whole-plant net production across the (prescribed midday) Tair gradient.
## TF24 (PM, no damage) vs TF24t with all thermal COSTS zeroed (isolates the
## jmax*N damage feedback + acclimation d_S shift) vs full TF24t (damage + costs).
zero_costs <- function(st) {
  st$c_acclim_maint <- 0; st$c_repair_maint <- 0; st$c_repair_flux <- 0
  st$c_build_topt <- 0; st$c_build_tcrit <- 0; st$c_accl_induct <- 0
  st
}
net_at <- function(type, Tmid, damage_only = FALSE) {
  env <- Environment(type)
  env$extrinsic_drivers_set_constant("air_temp", Tmid)
  if (type == "TF24") {
    st <- TF24_Strategy(); st$pars$use_energy_balance <- 1
    ind <- TF24_Individual(st)
  } else {
    st <- TF24t_Strategy(); if (damage_only) st <- zero_costs(st)
    ind <- TF24t_Individual(st)
  }
  ind$set_state("height", 8)
  ind$compute_rates(env)
  ind$net_mass_production_dt(env)
}
Tgrid <- c(25, 32, 38, 42, 46)
plant <- data.frame(
  Tair       = Tgrid,
  net_pm     = vapply(Tgrid, function(t) net_at("TF24", t), numeric(1)),
  net_dmg    = vapply(Tgrid, function(t) net_at("TF24t", t, damage_only = TRUE), numeric(1)),
  net_full   = vapply(Tgrid, function(t) net_at("TF24t", t), numeric(1)))
plant$rel_dmg_pct  <- 100 * (plant$net_dmg  - plant$net_pm) / abs(plant$net_pm)
plant$rel_full_pct <- 100 * (plant$net_full - plant$net_pm) / abs(plant$net_pm)

## ---------------------------------------------------------------------------
## (D) Acclimation buffering at a prescribed heatwave leaf temperature within
## acclimation's reach (Tleaf = 42; T_crit can rise to 38 + 6 = 44).
accl <- data.frame(A_crit = c(0, 0.5, 1, 2, 5))
cells <- lapply(accl$A_crit, function(a)
  solve_cell(1500, 42, 2, atls = TRUE, pm = FALSE, A_crit = a))
accl$N <- vapply(cells, function(x) x$N, numeric(1))
accl$A <- vapply(cells, function(x) x$A, numeric(1))

## ---------------------------------------------------------------------------
## Scorecard + gate decision.
cat("=== (A) CLEAN mechanism (prescribed Tleaf): damage OFF vs ON ===\n")
print(round(resA, 3), row.names = FALSE)
finA <- resA[is.finite(resA$relA_pct), ]
band <- subset(finA, Tleaf >= 38 & Tleaf <= 44)
cat(sprintf("max |relative dA| over Tleaf 38-44 = %.1f %%\n", max(abs(band$relA_pct))))
cat(sprintf("N at Tleaf=30 = %.3f ; Tleaf=40 = %.3f ; Tleaf=44 = %.3f\n",
            resA$N[resA$Tleaf==30][1], resA$N[resA$Tleaf==40][1], resA$N[resA$Tleaf==44][1]))

cat("\n=== (B) REALISTIC PM operating point ===\n")
print(round(resB, 3), row.names = FALSE)

cat("\n=== (C) whole-plant net production (height=8) ===\n")
print(round(plant, 4), row.names = FALSE)

cat("\n=== (D) T_crit acclimation buffers N/A at prescribed Tleaf=42 ===\n")
print(round(accl, 4), row.names = FALSE)

## Gate (mirror #523): the ATLS damage mechanism is MATERIAL iff the clean leaf
## response exceeds 5% in the damage band AND acclimation demonstrably buffers it.
gate_leaf  <- max(abs(band$relA_pct)) > 5
gate_accl  <- (accl$N[accl$A_crit == 5] - accl$N[accl$A_crit == 0]) > 0.05
gate_plant <- max(abs(plant$rel_dmg_pct[plant$Tair >= 38])) > 5
cat(sprintf("\nGATE: leaf-band |relA|>5%%: %s | acclim buffers N: %s | plant damage-only >5%%: %s\n",
            gate_leaf, gate_accl, gate_plant))
cat(sprintf("=> %s\n", if (gate_leaf && gate_accl) "KEEP (damage feedback material and acclimation-buffered)"
                       else "DROP (immaterial)"))
