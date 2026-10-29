## Doc 7.5.2 and 7.5.3: leaf-level sensitivity of the PM path to its two
## least-constrained inputs, the boundary-layer resistance ra and the net
## radiation Rn. Same leaf and same 64-cell factorial as step5_fick_vs_pm.R;
## every run has PM on, and each perturbation is injected into the leaf after
## set_physiology has derived ra/Rn and before the solve reads them.
##
## Every cell is reported at its GLOBAL profit maximum (pm_solve_cell with
## global = TRUE), not at the solver's answer: on the PM path the profit can
## have two maxima and the solver does not always return the higher, and a
## perturbation that merely moves the solver between them would read as a
## sensitivity of the model. The `gap` column records where that happened.
##
## Run from the plant repo root:
##   Rscript notes/penman-monteith/sensitivity_leaf.R [outdir]
## Writes sensitivity_ra.csv, sensitivity_ra_canopy.csv, sensitivity_rn.csv.
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
outdir <- if (length(args) >= 1) args[[1]] else here
source(file.path(here, "pm_leaf_common.R"))
pm_check_stack()
t0 <- proc.time()[["elapsed"]]

grid <- expand.grid(PAR = c(200, 500, 1000, 2000),
                    Tair = c(15, 25, 35, 40),
                    VPD = c(0.5, 1, 2, 3))
ra_default <- PM_C_RA * sqrt(0.05 / 2)   # wind model at d = 0.05 m, U0 = 2 m/s

## The injection is the whole method, so check it before using it: injecting
## the model's own ra changes nothing, and injecting ra = 50 is the same leaf as
## reaching ra = 50 through the wind model (U0 = 0.8 m/s).
chk <- expand.grid(PAR = c(500, 2000), Tair = c(25, 40), VPD = 2)
for (i in seq_len(nrow(chk))) {
  g <- chk[i, ]
  a <- pm_solve_cell(g$PAR, g$Tair, g$VPD, TRUE)
  b <- pm_solve_cell(g$PAR, g$Tair, g$VPD, TRUE, ra = a$ra_model)
  stopifnot(identical(a[c("Tleaf", "E", "A")], b[c("Tleaf", "E", "A")]))
  w <- pm_solve_cell(g$PAR, g$Tair, g$VPD, TRUE, U = 0.8)
  x <- pm_solve_cell(g$PAR, g$Tair, g$VPD, TRUE, ra = 50)
  stopifnot(abs(w$ra - 50) < 1e-12,
            isTRUE(all.equal(w[c("Tleaf", "E", "A")], x[c("Tleaf", "E", "A")],
                             tolerance = 1e-12)))
}
cat("injection checks passed\n")

## ---- 7.5.2a: fixed ra vs the wind model, across the factorial ---------------
## ra levels (s/m): 20 and 63 are the wind model at U0 = 5 and 0.5 m/s; 31.6 is
## the default (U0 = 2); 50 is the doc's and the leaf's fixed fallback; 110 is
## the doc's exponential profile one unit of LAI into the canopy (alpha_w = 2.5);
## 200 is the top of the doc's "realistic" 20-200 range.
ra_levels <- c(20, ra_default, 50, PM_C_RA * sqrt(0.05 / 0.5),
               ra_default * exp(2.5 * 1 / 2), 200)
ra_rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  for (ra in ra_levels) {
    ra_rows[[length(ra_rows) + 1L]] <-
      pm_solve_cell(g$PAR, g$Tair, g$VPD, TRUE, ra = ra, global = TRUE)
  }
}
ra_res <- do.call(rbind, ra_rows)
write.csv(ra_res, file.path(outdir, "sensitivity_ra.csv"), row.names = FALSE)

## ---- 7.5.2b: the within-canopy profile, along a depth transect --------------
## Cumulative LAI L from the crown top. Absorbed PAR follows TF24's own light
## model, k_I * PPFD * exp(-k_I * L) at the default PPFD = 1800 and k_I = 0.5.
## Three ra treatments: the wind model at the top value for every layer (what
## plant does now), the doc's profile ra0 * exp(alpha_w * L / 2) with
## alpha_w = 2.5, and the fixed 50 s/m.
k_I <- PM_PARS$k_I
transect <- expand.grid(L = c(0, 0.5, 1, 1.5, 2), Tair = c(15, 25, 35),
                        VPD = c(1, 2),
                        treatment = c("wind_top", "profile", "fixed50"),
                        stringsAsFactors = FALSE)
can_rows <- vector("list", nrow(transect))
for (i in seq_len(nrow(transect))) {
  g <- transect[i, ]
  PAR <- k_I * 1800 * exp(-k_I * g$L)
  ra <- switch(g$treatment, wind_top = ra_default,
               profile = ra_default * exp(2.5 * g$L / 2), fixed50 = 50)
  can_rows[[i]] <- cbind(L = g$L, treatment = g$treatment,
                         pm_solve_cell(PAR, g$Tair, g$VPD, TRUE, ra = ra, global = TRUE))
}
can_res <- do.call(rbind, can_rows)
write.csv(can_res, file.path(outdir, "sensitivity_ra_canopy.csv"), row.names = FALSE)

## ---- 7.5.3: net radiation +-30 %, and the longwave term on its own -----------
## (a) Rn scaled by 0.7 and 1.3 (the doc's +-30 %), whole Rn including the
## fixed -40 W/m2 longwave offset.
## (b) The longwave term alone: replace the fixed -40 W/m2 with the doc's
## clear-sky form eps_l * sigma * Ta^4 * (eps_a - 1), eps_l = 0.97, at
## eps_a = 0.7, 0.8, 0.9 (doc 3.3 and 6), keeping the shortwave 2*PAR/4.57.
sigma <- 5.67e-8
lw_net <- function(Tair, eps_a) 0.97 * sigma * (Tair + 273.15)^4 * (eps_a - 1)
rn_rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  for (s in c(0.7, 1.0, 1.3)) {
    rn_rows[[length(rn_rows) + 1L]] <- cbind(
      variant = sprintf("Rn_x%.1f", s), LW = NA_real_,
      pm_solve_cell(g$PAR, g$Tair, g$VPD, TRUE, rn_scale = s, global = TRUE))
  }
  sw <- 2 * g$PAR / 4.57
  for (ea in c(0.7, 0.8, 0.9)) {
    lw <- lw_net(g$Tair, ea)
    base_rn <- sw - 40
    rn_rows[[length(rn_rows) + 1L]] <- cbind(
      variant = sprintf("LW_eps%.1f", ea), LW = lw,
      pm_solve_cell(g$PAR, g$Tair, g$VPD, TRUE,
                    rn_scale = (sw + lw) / base_rn, global = TRUE))
  }
}
rn_res <- do.call(rbind, rn_rows)
write.csv(rn_res, file.path(outdir, "sensitivity_rn.csv"), row.names = FALSE)

cat(sprintf("all leaf solves done in %.1f s\n", proc.time()[["elapsed"]] - t0))

## ---- Summaries -----------------------------------------------------------------
## A cell is SHUT when its base PM operating point moves no water (E < 1e-6):
## the leaf sits at its no-transpiration temperature and A = -R_d.
## Relative A is reported only where the base A > 1 umol m-2 s-1.
contrast <- function(base, pert, label) {
  shut <- base$E < 1e-6
  relA <- ifelse(base$A > 1, 100 * (pert$A - base$A) / base$A, NA_real_)
  relE <- ifelse(base$E > 1e-6, 100 * (pert$E - base$E) / base$E, NA_real_)
  dT <- pert$Tleaf - base$Tleaf
  cat(sprintf(paste0("%-24s |dTleaf| median %.2f max %.2f C; ",
                     "|dA| median %.1f max %.1f %% (n=%d); ",
                     "|dE| median %.1f max %.1f %%; cells changing shut state: %d; ",
                     "solver missed the global max (>0.1) in %d\n"),
              label, median(abs(dT)), max(abs(dT)),
              median(abs(relA), na.rm = TRUE), max(abs(relA), na.rm = TRUE),
              sum(is.finite(relA)),
              median(abs(relE), na.rm = TRUE), max(abs(relE), na.rm = TRUE),
              sum(shut != (pert$E < 1e-6)), sum(pert$gap > 0.1)))
}
key <- function(d) paste(d$PAR, d$Tair, d$VPD)
pick <- function(d, keep) { d <- d[keep, ]; d[order(key(d)), ] }

cat("--- 7.5.2a: ra = 50 (fixed) against the wind model's 31.6, PM on ---\n")
base <- pick(ra_res, abs(ra_res$ra - ra_default) < 1e-9)
for (r in ra_levels[abs(ra_levels - ra_default) > 1e-9]) {
  contrast(base, pick(ra_res, abs(ra_res$ra - r) < 1e-9), sprintf("ra %.1f vs 31.6", r))
}
cat(sprintf("ra produced by the wind model over U0 in [0.5, 5] m/s: %.1f-%.1f s/m\n",
            PM_C_RA * sqrt(0.05 / 5), PM_C_RA * sqrt(0.05 / 0.5)))
cat(sprintf("ra produced by the doc profile (alpha_w 2.5-3.5) at L = 1, 2, 3: %s s/m\n",
            paste(sprintf("%.0f-%.0f", ra_default * exp(2.5 * 1:3 / 2),
                          ra_default * exp(3.5 * 1:3 / 2)), collapse = ", ")))

cat("--- 7.5.2b: canopy transect (Tleaf, A by depth) ---\n")
print(reshape(can_res[, c("L", "Tair", "VPD", "treatment", "Tleaf", "A")],
              idvar = c("L", "Tair", "VPD"), timevar = "treatment",
              direction = "wide"), row.names = FALSE, digits = 3)
## Crown-integrated A over L in [0, 2], trapezium, profile vs wind_top.
crown <- aggregate(A ~ Tair + VPD + treatment, can_res, function(a) {
  L <- c(0, 0.5, 1, 1.5, 2); sum(diff(L) * (head(a, -1) + tail(a, -1)) / 2)
})
print(reshape(crown, idvar = c("Tair", "VPD"), timevar = "treatment",
              direction = "wide"), row.names = FALSE, digits = 4)

cat("--- 7.5.3: Rn +-30 % and the longwave term, PM on ---\n")
base <- pick(rn_res, rn_res$variant == "Rn_x1.0")
for (v in c("Rn_x0.7", "Rn_x1.3", "LW_eps0.7", "LW_eps0.8", "LW_eps0.9")) {
  contrast(base, pick(rn_res, rn_res$variant == v), v)
}
lw <- outer(c(15, 25, 35, 40), c(0.7, 0.8, 0.9), lw_net)
dimnames(lw) <- list(Tair = c(15, 25, 35, 40), eps_a = c(0.7, 0.8, 0.9))
cat("LW_net (W/m2) by Tair x eps_a, against the fixed -40:\n"); print(round(lw, 1))
