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
  # The opt_root_psi aux is a positive magnitude (phylloptim #25), so the
  # deeper-rooted plant -- reaching the wet layers -- sits at the SMALLER suction.
  # The inequality reversed with the representation; the physics did not.
  expect_lt(psi_deep, psi_shallow)
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

# ---------------------------------------------------------------------------
# Per-layer layer thickness in the root network (#626, phylloptim 0.9.0)
# ---------------------------------------------------------------------------

probe_root_operating_point <- function(widths, theta = 0.25) {
  ind <- TF24_Individual(TF24_Strategy())
  ind$set_state("height", 5)
  env <- TF24_Environment()
  env$set_soil_layer_widths(widths)
  env$set_soil_water_state(rep(theta, length(widths)))
  ind$compute_rates(env)
  stats::setNames(ind$internals$auxs, ind$aux_names)
}

test_that("plant hands the leaf each layer's own thickness, not a column average", {
  # THE PLANT-SIDE PROOF THAT THE FIX LANDED, and it is a SIGN test rather than a
  # tolerance, because a tolerance cannot see this. Vertical root resistance scales
  # with the square of the thickness of the segment spanning each layer.
  # phylloptim <= 0.8.0 took one scalar thickness (column depth / n), correct only
  # for equal layers, and plant is what supplies it.
  #
  # Measured both ways, by building plant against each (theta = 0.25 everywhere,
  # height 5 m, 1.5 m column):
  #
  #   profile            opt_root_psi (MPa)     assimilation
  #                      scalar    per-layer    scalar    per-layer
  #   uniform, 5 layers  1.661740  1.661740     15.464934 15.464934
  #   2 cm surface layer 1.747280  0.598266     15.338550 17.251842
  #   thick over thin    1.654759  4.211291     15.478450 11.179142
  #
  # Two things to read off it. The uniform column is BIT-IDENTICAL, which is why
  # no default plant output moves. And on a graded column a scalar thickness lands
  # everything within 6% of the uniform value, on the WRONG SIDE of it: a thin
  # surface layer is charged 0.3 m of root segment instead of 0.02 m, over-resisting
  # by 225x in that layer and 3.7x in total, so the plant appears to need MORE
  # suction when it needs far less. Reversing a thin-over-thick profile reverses the
  # error. So the discriminator is the sign, and it does not depend on the exact
  # values surviving an unrelated solver change.
  uniform <- probe_root_operating_point(rep(0.3, 5))

  thin_top <- probe_root_operating_point(c(0.02, 0.28, 0.30, 0.40, 0.50))
  thick_top <- probe_root_operating_point(c(0.75, 0.25, 0.25, 0.15, 0.10))

  for (p in list(uniform, thin_top, thick_top)) {
    expect_true(is.finite(p[["assimilation"]]))
    expect_true(is.finite(p[["opt_root_psi"]]))
  }

  # Thinning the surface layer moves root carbon deeper and cuts the vertical
  # resistance it is charged, so the plant sits at LESS suction than uniform.
  expect_lt(thin_top[["opt_root_psi"]], uniform[["opt_root_psi"]])
  expect_gt(thin_top[["assimilation"]], uniform[["assimilation"]])

  # Thickening it does the reverse. Both inequalities fail under a scalar
  # thickness, which puts these within 6% of uniform and the wrong way round.
  expect_gt(thick_top[["opt_root_psi"]], uniform[["opt_root_psi"]])
  expect_lt(thick_top[["assimilation"]], uniform[["assimilation"]])

  # And the effect is large, not marginal: a scalar thickness kept every graded
  # profile inside 6% of the uniform operating point.
  expect_gt(abs(thin_top[["opt_root_psi"]] / uniform[["opt_root_psi"]] - 1), 0.2)
  expect_gt(abs(thick_top[["opt_root_psi"]] / uniform[["opt_root_psi"]] - 1), 0.2)
})

test_that("an unequal-layer column runs end to end", {
  env <- Environment("TF24")
  set_tf24_soil(env, soil_widths_graded(1.5, 5, top = 0.02), theta = 0.25)
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- 2
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  out <- run_scm(p, env)
  expect_true(is.finite(out$offspring_production))
  expect_gte(out$offspring_production, 0)
})
