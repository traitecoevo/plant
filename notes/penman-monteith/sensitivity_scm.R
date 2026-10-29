## Doc 7.5.4: stand-level runs, PM on vs off, at several air temperatures.
## TF24's SCM with constant drivers at their defaults except `leaf_temp`, which
## is the air temperature here (the PM path reinterprets it as such). The
## density is carried in birth date (Control$node_density_in_birth_date =
## TRUE), the coordinate TF24 needs (#590).
##
## Up to three runs per scenario:
##   y1  one species, a 1-year patch: what a year of growth builds (census at
##       t = 1: leaf area, above-ground mass) and the water it takes up;
##   p5  the same species over a 5-year patch: offspring production;
##   c5  two species (lma 0.0825 vs 0.10) over a 5-year patch: offspring of
##       each, i.e. whether PM changes who wins (main grid only).
## Traits and birth rates are those of test-strategy-tf24.R's 5-year fixtures.
##
## Scenarios:
##   main     Tair in {10..35}, PM off and on, every driver else at default;
##   matched  Fick at PM's canopy-top leaf temperature (from the main y1 run),
##            once at the air deficit (1 kPa) and once at the leaf-to-air
##            deficit PM's top leaf sees -- separates "the leaf is hotter"
##            from "the leaf sees a bigger deficit";
##   wind     PM on at U0 = 0.8 and 5 m/s, i.e. ra = 50 (the fixed fallback)
##            and 20 s/m, against the default U0 = 2 (ra = 31.6): 7.5.2 at
##            stand level.
##
## Run from the plant repo root:
##   Rscript notes/penman-monteith/sensitivity_scm.R [outdir]
## Writes sensitivity_scm.csv. ~12 min unloaded, 10 of them in the one Tair = 35
## PM-on cell, whose stand dies (that cell is CPU-bound: 61 min under heavy load).
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
outdir <- if (length(args) >= 1) args[[1]] else here
source(file.path(here, "pm_leaf_common.R"))
pm_check_stack()

ctrl <- Control(node_density_in_birth_date = TRUE)
stopifnot(isTRUE(ctrl$node_density_in_birth_date))

## The switch lives in each strategy's pars. Assigning it anywhere else (for
## example at the strategy's top level) is silently ignored and reproduces the
## default arm, so it is asserted after it is set, on every strategy.
set_pm <- function(p, on) {
  for (i in seq_along(p$strategies)) {
    s <- p$strategies[[i]]
    s$pars$use_energy_balance <- on
    p$strategies[[i]] <- s
  }
  got <- vapply(p$strategies, function(s) s$pars$use_energy_balance, 0)
  stopifnot(all(got == on))
  p
}

make_p <- function(lma, lifetime) {
  p0 <- scm_base_parameters("TF24")
  p0$max_patch_lifetime <- lifetime
  add_strategies(p0, trait_matrix(c(lma, rep(5, length(lma))),
                                  c("lma", "hmat")),
                 hyperpar = TF24_hyperpar,
                 birth_rate = as.list(rep(20, length(lma))))
}

env0 <- Environment("TF24")
drv <- vapply(env0$extrinsic_drivers_get_names(),
              function(n) env0$extrinsic_drivers_evaluate(n, 0), 0)
cat("driver defaults:", paste(sprintf("%s=%g", names(drv), drv), collapse = ", "), "\n")

make_env <- function(leaf_temp, vpd = drv[["atm_vpd"]], wind = drv[["wind_speed"]]) {
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("leaf_temp", leaf_temp)
  env$extrinsic_drivers_set_constant("atm_vpd", vpd)
  env$extrinsic_drivers_set_constant("wind_speed", wind)
  env
}

run_one <- function(p, env, pm) {
  t0 <- proc.time()[["elapsed"]]
  scm <- run_scm(set_pm(p, pm), env, ctrl, run = FALSE)
  scm$collect <- TRUE
  scm$run()
  res <- plant:::tidy_patch(lapply(scm$history, "[[", "state"))
  sp <- res$species
  last <- sp[sp$time == max(sp$time) & is.finite(sp$height), ]
  top <- last[which.max(last$height), ]
  flux <- res$env$soil_moist_cumulative_flux
  flux <- flux[nrow(flux), ]
  cen <- stand_census(scm)
  list(elapsed = proc.time()[["elapsed"]] - t0,
       offspring = scm$offspring_production,
       leaf_area = cen[["leaf_area"]], mass_above_ground = cen[["mass_above_ground"]],
       water_uptake = flux$sum_resource_depletion,
       top_height = top$height, top_Tleaf = top$Tleaf, top_A = top$assimilation,
       top_E = top$transpiration, top_psi_stem = top$opt_psi_stem,
       top_net_dt = top$net_mass_production_dt)
}

## One scenario: air temperature `Tair` (for the record), the driver values
## actually set, PM on or off, and whether to run the two-species patch.
run_scenario <- function(scenario, Tair, pm, leaf_temp = Tair,
                         vpd = drv[["atm_vpd"]], wind = drv[["wind_speed"]],
                         comp = FALSE) {
  env <- function() make_env(leaf_temp, vpd, wind)
  y1 <- run_one(make_p(0.0825, 1), env(), pm)
  p5 <- run_one(make_p(0.0825, 5), env(), pm)
  c5 <- if (comp) run_one(make_p(c(0.0825, 0.10), 5), env(), pm)
        else list(offspring = c(NA_real_, NA_real_), elapsed = 0)
  r <- data.frame(scenario = scenario, Tair = Tair, pm = pm,
                  leaf_temp_driver = leaf_temp, atm_vpd = vpd, wind_speed = wind,
                  y1_leaf_area = y1$leaf_area, y1_mass_above_ground = y1$mass_above_ground,
                  y1_water_uptake = y1$water_uptake,
                  y1_top_height = y1$top_height, y1_top_Tleaf = y1$top_Tleaf,
                  y1_top_A = y1$top_A, y1_top_E = y1$top_E,
                  y1_top_psi_stem = y1$top_psi_stem, y1_top_net_dt = y1$top_net_dt,
                  p5_offspring = p5$offspring, p5_leaf_area = p5$leaf_area,
                  p5_mass_above_ground = p5$mass_above_ground,
                  p5_water_uptake = p5$water_uptake, p5_top_height = p5$top_height,
                  p5_top_Tleaf = p5$top_Tleaf,
                  c5_offspring_lma0825 = c5$offspring[[1]],
                  c5_offspring_lma10 = c5$offspring[[2]],
                  seconds = y1$elapsed + p5$elapsed + c5$elapsed)
  cat(sprintf(paste0("%-8s Tair %2g pm %d leaf_temp %5.2f vpd %4.2f wind %3.1f: ",
                     "y1 mass %.4g leafA %.4g water %.4g topTleaf %.2f | ",
                     "p5 offspring %.4g | c5 %.4g / %.4g | %.0f s\n"),
              scenario, Tair, pm, leaf_temp, vpd, wind, r$y1_mass_above_ground,
              r$y1_leaf_area, r$y1_water_uptake, r$y1_top_Tleaf, r$p5_offspring,
              r$c5_offspring_lma0825, r$c5_offspring_lma10, r$seconds))
  r
}

rows <- list()
for (Tair in c(10, 15, 20, 25, 30, 35)) for (pm in c(0, 1)) {
  rows[[length(rows) + 1L]] <- run_scenario("main", Tair, pm, comp = TRUE)
}
main <- do.call(rbind, rows)

## Matched: only where PM's stand survives its first year.
for (Tair in c(10, 15, 20, 25)) {
  Tl <- main$y1_top_Tleaf[main$scenario == "main" & main$Tair == Tair & main$pm == 1]
  vpd_leaf <- drv[["atm_vpd"]] + pm_esat(Tl) - pm_esat(Tair)
  rows[[length(rows) + 1L]] <- run_scenario("matchT", Tair, 0, leaf_temp = Tl)
  rows[[length(rows) + 1L]] <- run_scenario("matchTV", Tair, 0, leaf_temp = Tl,
                                            vpd = vpd_leaf)
}
for (Tair in c(10, 15, 20, 25)) for (U in c(0.8, 5)) {
  rows[[length(rows) + 1L]] <- run_scenario("wind", Tair, 1, wind = U)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(outdir, "sensitivity_scm.csv"), row.names = FALSE)

## PM relative to Fick, per air temperature (main grid).
off <- main[main$pm == 0, ]; on <- main[main$pm == 1, ]
rel <- function(a, b) 100 * (b - a) / a
share <- function(d) d$c5_offspring_lma0825 / (d$c5_offspring_lma0825 + d$c5_offspring_lma10)
cmp <- data.frame(Tair = off$Tair,
                  dTleaf_top_y1 = on$y1_top_Tleaf - off$Tair,
                  mass_y1_pct = rel(off$y1_mass_above_ground, on$y1_mass_above_ground),
                  leaf_area_y1_pct = rel(off$y1_leaf_area, on$y1_leaf_area),
                  water_y1_pct = rel(off$y1_water_uptake, on$y1_water_uptake),
                  mass_p5_pct = rel(off$p5_mass_above_ground, on$p5_mass_above_ground),
                  offspring_ratio_p5 = on$p5_offspring / off$p5_offspring,
                  comp_share_lma0825_off = share(off), comp_share_lma0825_on = share(on))
cat("--- main: PM relative to Fick ---\n")
print(signif(cmp, 4), row.names = FALSE)
cat("--- matched and wind scenarios ---\n")
print(res[res$scenario != "main",
          c("scenario", "Tair", "pm", "leaf_temp_driver", "atm_vpd", "wind_speed",
            "y1_mass_above_ground", "y1_water_uptake", "y1_top_Tleaf",
            "p5_offspring")], digits = 4, row.names = FALSE)
cat(sprintf("total runtime %.0f s\n", sum(res$seconds)))
