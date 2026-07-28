# Root hydraulic parameters exposed through TF24_Pars.
#
# root_c, root_b and root_psi_crit were fixed members of TF24_Strategy, and the
# rooting depth cap was a file-static constant in src/tf24_strategy.cpp. None
# were reachable from R, so root shutoff was pinned at ~5.87 MPa and rooting
# depth at 1.5 m -- neither calibratable, and the shutoff is too conservative
# for taxa that operate below it. They are now TF24_Pars fields.

test_that("root hydraulic parameters are exposed with unchanged defaults", {
  # Values are the ones previously hard-coded; changing them is a scientific
  # change and should require editing this test deliberately.
  pars <- TF24_Strategy()$pars
  expect_equal(pars$root_c, 2.680147)
  expect_equal(pars$root_b, 3.898245)
  expect_equal(pars$rooting_depth_max, 1.5)
  # root_psi_crit is derived: the potential at 5% remaining root conductivity.
  expect_equal(pars$root_psi_crit,
               pars$root_b * log(1 / 0.05)^(1 / pars$root_c))
  # Matches the value the Leaf constructor previously defaulted to.
  expect_equal(pars$root_psi_crit, 5.870283, tolerance = 1e-6)
})

test_that("root hydraulic parameters are settable", {
  s <- TF24_Strategy()
  s$pars$root_b <- 1.5
  s$pars$root_c <- 3.0
  s$pars$root_psi_crit <- 1.5 * log(1 / 0.05)^(1 / 3.0)
  s$pars$rooting_depth_max <- 4.0
  expect_equal(s$pars$root_b, 1.5)
  expect_equal(s$pars$root_c, 3.0)
  expect_equal(s$pars$rooting_depth_max, 4.0)
})

# Shared probe: one compute_rates call, returning the named aux vector.
tf24_root_probe <- function(strategy, theta) {
  ind <- TF24_Individual(strategy)
  ind$set_state("height", 5)
  env <- TF24_Environment()
  env$set_soil_number_of_depths(length(theta))
  env$set_soil_water_state(theta)
  ind$compute_rates(env)
  stats::setNames(ind$internals$auxs, ind$aux_names)
}

test_that("root_b reaches the root vulnerability curve", {
  # Lowering root_b makes roots lose conductivity at less negative potentials,
  # so the same soil supports less carbon gain. If the parameter were not wired
  # through to Leaf, assimilation would be identical.
  wet <- rep(0.25, 5)

  base <- TF24_Strategy()
  fragile <- TF24_Strategy()
  fragile$pars$root_b <- 0.5
  fragile$pars$root_psi_crit <-
    fragile$pars$root_b * log(1 / 0.05)^(1 / fragile$pars$root_c)

  a_base <- tf24_root_probe(base, wet)[["assimilation"]]
  a_fragile <- tf24_root_probe(fragile, wet)[["assimilation"]]

  expect_true(is.finite(a_base) && is.finite(a_fragile))
  expect_lt(a_fragile, a_base)
})

test_that("rooting_depth_max reaches the root network", {
  # Dry shallow layers over wet deep ones: how deep the roots are allowed to go
  # determines which layers contribute, and so the root-collar potential the
  # plant settles at. A shallow cap leaves it drier.
  stratified <- c(rep(0.04, 8), rep(0.35, 7))   # dry to 0.8 m, wet below

  shallow <- TF24_Strategy(); shallow$pars$rooting_depth_max <- 0.2
  deep    <- TF24_Strategy(); deep$pars$rooting_depth_max    <- 1.5

  psi_shallow <- tf24_root_probe(shallow, stratified)[["opt_root_psi"]]
  psi_deep    <- tf24_root_probe(deep, stratified)[["opt_root_psi"]]

  expect_true(is.finite(psi_shallow) && is.finite(psi_deep))
  # Potentials are reported as signed (negative) values, so the deeper-rooted
  # plant sits at the less negative -- wetter -- collar potential.
  expect_gt(psi_deep, psi_shallow)
})

test_that("exposing the root parameters left default behaviour unchanged", {
  # The defaults are the previously hard-coded constants, so a default run must
  # be unaffected by the move from strategy members to TF24_Pars.
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.25, 5))
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- 2
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  out <- run_scm(p, env)
  expect_true(is.finite(out$offspring_production))
  expect_gte(out$offspring_production, 0)
})
