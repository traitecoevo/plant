# Smoke test for the staging Penman-Monteith leaf demo (#523). Exercises the
# exact leaf-driving helpers the overstorey_staging/ vignette uses, so the demo
# cannot silently rot. overstorey_staging/ is .Rbuildignore'd (not installed), so
# this is dev-only: skip when the helper file is absent (installed/CRAN checks).

test_that("PM leaf demo helpers run end-to-end for Fick and PM", {
  skip_on_cran()
  helpers <- test_path("..", "..", "overstorey_staging", "pm_demo_helpers.R")
  skip_if_not(file.exists(helpers), "overstorey_staging/ not present (built package)")
  source(helpers, local = TRUE)

  # A small environment grid spanning cool -> hot / high-radiation.
  grid <- expand.grid(PAR = c(400, 2000), Tair = c(20, 40), VPD = c(1, 3))
  res <- pm_solve_grid(grid)

  # Both models solved for every cell, all outputs finite.
  expect_equal(nrow(res), nrow(grid) * 2L)
  expect_true(all(c(FALSE, TRUE) %in% res$pm))
  for (v in c("Tleaf", "opt_psi_stem", "A", "gs", "E", "profit")) {
    expect_true(all(is.finite(res[[v]])),
                info = paste("non-finite", v, "in demo grid"))
  }

  # Fick keeps Tleaf == Tair; PM departs from it under high radiation.
  fick <- subset(res, !pm)
  expect_equal(fick$Tleaf, fick$Tair)
  pm <- subset(res, pm)
  hot <- subset(pm, PAR == 2000 & Tair == 40)
  expect_true(all(hot$Tleaf > hot$Tair))          # net warming when hot + bright

  # Profit-anatomy scan over root-collar potential (the actual solver decision
  # variable, see pm_collar_curve) stays finite for both models.
  for (use_pm in c(FALSE, TRUE)) {
    curve <- pm_collar_curve(1800, 38, 2.5, use_pm, seq(0.31, 5.5, length.out = 20))
    expect_true(all(is.finite(curve$profit)))
    expect_true(all(is.finite(curve$assim)))
  }
})

test_that("PM path fails fast on a non-finite wind speed (review: itowers1)", {
  skip_on_cran()
  helpers <- test_path("..", "..", "overstorey_staging", "pm_demo_helpers.R")
  skip_if_not(file.exists(helpers), "overstorey_staging/ not present (built package)")
  source(helpers, local = TRUE)

  set_phys <- function(l) {
    l$set_physiology( root_carbon_per_leaf_area = 20,
                     PPFD = 1000, psi_soil = 0.3, soil_depth = 1,
                     leaf_specific_conductance_max = 5e-3, atm_vpd = 2, ca = 40, leaf_temp = 30,
                     atm_o2_kpa = 21, atm_kpa = 101.3)
  }

  # PM on + NA wind speed is a broken driver -> fail fast (not silent fallback).
  l <- pm_make_leaf(); l$use_energy_balance_ <- TRUE; l$d_ <- 0.05
  l$wind_speed_ <- NA_real_
  expect_error(set_phys(l), "non-finite wind_speed")

  # PM off + NA wind speed is fine: the wind model is not read off the PM path.
  l2 <- pm_make_leaf(); l2$use_energy_balance_ <- FALSE; l2$wind_speed_ <- NA_real_
  expect_no_error(set_phys(l2))

  # Zero wind (physically ra -> infinity) is a legitimate case, not an error:
  # it falls back to the fixed ra.
  l3 <- pm_make_leaf(); l3$use_energy_balance_ <- TRUE; l3$wind_speed_ <- 0
  expect_no_error(set_phys(l3))
})
