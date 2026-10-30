## Step 5 (doc 7.5.1): factorial Fick vs Penman-Monteith at the leaf level.
## Drives a single TF24-default leaf across Tair x VPD x PAR, PM off vs on, and
## reports dTleaf, dE, dA, dgs. Writes a scorecard RDS + CSV.
##
## Run from the plant repo root:
##   Rscript notes/penman-monteith/step5_fick_vs_pm.R <outstem>
##
## 2026-10 port: the design is July's (same grid, same leaf geometry, same
## outputs); only the API moved. The leaf lives in phylloptim now (#591), so the
## Leaf constructor and set_physiology take phylloptim's arguments, and the leaf
## is seated on TF24_floor at lambda_o = 0 exactly as TF24_Strategy seats it.
## The shared driver is pm_leaf_common.R.
args <- commandArgs(trailingOnly = TRUE)
outstem <- if (length(args) >= 1) args[[1]] else "step5"
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "pm_leaf_common.R"))
pm_check_stack()

solve_cell <- function(PAR, Tair, VPD, pm) pm_solve_cell(PAR, Tair, VPD, pm)

grid <- expand.grid(PAR = c(200, 500, 1000, 2000),
                    Tair = c(15, 25, 35, 40),
                    VPD = c(0.5, 1, 2, 3))
rows <- vector("list", nrow(grid))
t0 <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  off <- solve_cell(g$PAR, g$Tair, g$VPD, FALSE)
  on  <- solve_cell(g$PAR, g$Tair, g$VPD, TRUE)
  # Added in the 2026-10 port: the same cell at the GLOBAL profit maximum (a
  # brute-force scan of the collar bracket, pm_scan_collar). The columns above
  # stay the solver's, because that is what TF24 runs on; these say where it
  # did not find the higher of two maxima.
  off_g <- pm_solve_cell(g$PAR, g$Tair, g$VPD, FALSE, global = TRUE)
  on_g  <- pm_solve_cell(g$PAR, g$Tair, g$VPD, TRUE, global = TRUE)
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
    profit_off = off$profit, profit_on = on$profit,
    # Added in the 2026-10 port, to explain what moved since July: the
    # leaf-to-air deficit the PM leaf actually sees (#93) and dark respiration
    # at each path's operating point (#86).
    vpd_leaf_pm = on$vpd_leaf, Rd_off = off$R_d, Rd_on = on$R_d,
    gap_off = off_g$gap, gap_on = on_g$gap,
    A_off_global = off_g$A, A_on_global = on_g$A,
    E_on_global = on_g$E, Tleaf_pm_global = on_g$Tleaf)
}
res <- do.call(rbind, rows)
cat(sprintf("solved %d cells x 2 in %.1f s\n", nrow(res), proc.time()[["elapsed"]] - t0))

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
cat("--- cells where the solver missed the global profit maximum by > 0.1 ---\n")
miss <- subset(res, gap_off > 0.1 | gap_on > 0.1)
print(round(miss[, c("PAR","Tair","VPD","gap_off","gap_on","A_off","A_off_global",
                     "A_on","A_on_global","Tleaf_pm","Tleaf_pm_global")], 3),
      row.names = FALSE)
relA_g <- ifelse(res$A_off_global > 1,
                 100 * (res$A_on_global - res$A_off_global) / res$A_off_global, NA)
cat(sprintf("at the global maxima: max |Tleaf-Tair| = %.2f C, median %.2f C, max |relA| = %.1f %%\n",
            max(abs(res$Tleaf_pm_global - res$Tair)),
            median(abs(res$Tleaf_pm_global - res$Tair)), max(abs(relA_g), na.rm = TRUE)))
