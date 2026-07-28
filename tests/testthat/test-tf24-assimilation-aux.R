# The `assimilation` auxiliary variable.
#
# `assimilation` was listed in TF24_Strategy::aux_names() but never written: no
# aux_idx_assimilation existed and no set_aux call referenced it, so the slot
# reported whatever Internals happened to hold. That made carbon uptake
# unavailable as a model output. These tests pin the slot's content and units.
#
# Units: umol CO2 m^-2 leaf s^-1, and NET of dark respiration, because
# Leaf::assim_colimited() subtracts R_d_. Gross = assimilation + R_d_, with
# R_d_ = 0.015 * vcmax_ at the acclimated vcmax_.

tf24_aux <- function(strategy = TF24_Strategy(), height = 5, ppfd = 1800,
                     theta = rep(0.25, 5)) {
  ind <- TF24_Individual(strategy)
  ind$set_state("height", height)
  env <- TF24_Environment()
  env$extrinsic_drivers_set_constant("PPFD", ppfd)
  env$set_soil_number_of_depths(length(theta))
  env$set_soil_water_state(theta)
  ind$compute_rates(env)
  stats::setNames(ind$internals$auxs, ind$aux_names)
}

test_that("assimilation aux is written, not left unset", {
  aux <- tf24_aux()
  expect_true("assimilation" %in% names(aux))
  expect_true(is.finite(aux[["assimilation"]]))
  # A well-lit, well-watered plant must be fixing carbon. The pre-fix slot was
  # never assigned, so this is the assertion that would have caught it.
  expect_gt(aux[["assimilation"]], 0)
})

test_that("assimilation exceeds profit by the hydraulic cost", {
  # profit = assimilation - hydraulic_cost_TF(psi_stem), and the cost is
  # strictly positive whenever the plant is transpiring, so this ordering is
  # the invariant tying the two slots together.
  aux <- tf24_aux()
  expect_gt(aux[["assimilation"]], aux[["profit"]])
})

test_that("assimilation increases with light", {
  by_light <- vapply(c(0, 200, 800, 1800),
                     function(q) tf24_aux(ppfd = q)[["assimilation"]],
                     numeric(1))
  expect_true(all(diff(by_light) > 0))
})

test_that("assimilation in the dark is minus dark respiration", {
  # With no light there is no gross assimilation, so the reported net rate is
  # exactly -R_d_ = -0.015 * vcmax_. At the default air_temp of 25 C the
  # acclimated vcmax_ equals vcmax_25, so this is an exact expectation and
  # doubles as a check that the units are per unit leaf area, not canopy total.
  s <- TF24_Strategy()
  expected <- -0.015 * s$pars$vcmax_25
  expect_equal(tf24_aux(strategy = s, ppfd = 0)[["assimilation"]], expected,
               tolerance = 1e-8)
})

test_that("assimilation is finite when the hydraulic solve shuts down", {
  # Both shut-down exits in Leaf::find_root_collar_psi bypass
  # profit_psi_stem_TF, so assim_colimited_ used to be left at whatever the last
  # probe wrote. They now set it explicitly, keeping
  # profit == assimilation - hydraulic cost true in every branch.
  # Dry top layers over wet deep ones drives the shutdown.
  aux <- tf24_aux(theta = c(rep(0.04, 8), rep(0.35, 7)))
  expect_true(is.finite(aux[["assimilation"]]))
  expect_lt(aux[["assimilation"]], 0)
  s <- TF24_Strategy()
  expect_equal(aux[["assimilation"]], -0.015 * s$pars$vcmax_25, tolerance = 1e-8)
})

test_that("assimilation is reported under every shading model", {
  # The deep-crown model integrates each leaf output to a leaf-area-weighted
  # crown mean; assimilation must be integrated alongside the others rather than
  # left at the last quadrature point's value.
  for (model in c("mean-light", "crown-centre", "deep-crown")) {
    s <- TF24_Strategy()
    s$control$shading_model <- model
    aux <- tf24_aux(strategy = s)
    expect_true(is.finite(aux[["assimilation"]]), info = model)
    expect_gt(aux[["assimilation"]], 0)
    expect_gt(aux[["assimilation"]], aux[["profit"]])
  }
})
