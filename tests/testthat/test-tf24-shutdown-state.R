# Leaf state on the hydraulic shut-down path.
#
# `Leaf` is a value member of TF24_Model, reused across every compute_rates
# call for an individual. Both shut-down exits in find_root_collar_psi set
# profit_ directly and bypass profit_psi_stem_TF, and the first of them returns
# before any E_from_Soil_to_Root_Collar call in that solve. So every leaf output
# they do not assign keeps the *previous* step's value.
#
# That was not merely a reporting problem: soil_consumption_ feeds
# TF24_Model::evapotranspiration_dt and hence the patch water balance, so a
# plant whose stomata had closed carried on drawing its last wet-step uptake out
# of the soil. The water budget still closes in that state (what is recorded as
# depleted is what is removed), so the conservation tests cannot catch it --
# only these can.

tf24_aux_at <- function(individual, env, theta) {
  env$set_soil_water_state(rep(theta, env$get_soil_number_of_depths()))
  individual$compute_rates(env)
  stats::setNames(individual$internals$auxs, individual$aux_names)
}

test_that("shut-down zeroes the whole transport chain", {
  ind <- TF24_Individual(TF24_Model())
  ind$set_state("height", 5)
  env <- TF24_Environment()
  env$set_soil_number_of_depths(5)

  # A wet step first, so there is a non-zero state available to go stale.
  wet <- tf24_aux_at(ind, env, 0.30)
  expect_gt(wet[["transpiration"]], 0)
  expect_gt(wet[["E_up_"]], 0)
  expect_gt(wet[["stom_cond_CO2"]], 0)

  # Now drier than psi_crit, which triggers the first shut-down exit.
  dry <- tf24_aux_at(ind, env, 0.02)
  expect_equal(dry[["transpiration"]], 0)
  expect_equal(dry[["stom_cond_CO2"]], 0)
  expect_equal(dry[["E_up_"]], 0)
  # The regression: these used to equal the wet-step values exactly.
  expect_false(isTRUE(all.equal(dry[["transpiration"]], wet[["transpiration"]])))
  expect_false(isTRUE(all.equal(dry[["E_up_"]], wet[["E_up_"]])))
})

test_that("shut-down is not a one-way door", {
  # Zeroing must not corrupt the leaf state: rewetting has to recover the same
  # operating point as an equivalent plant that never dried.
  ind <- TF24_Individual(TF24_Model())
  ind$set_state("height", 5)
  env <- TF24_Environment()
  env$set_soil_number_of_depths(5)

  before <- tf24_aux_at(ind, env, 0.30)
  invisible(tf24_aux_at(ind, env, 0.02))   # shut down
  after <- tf24_aux_at(ind, env, 0.30)     # rewet

  for (v in c("transpiration", "E_up_", "stom_cond_CO2", "assimilation",
              "profit", "opt_psi_leaf")) {
    expect_equal(after[[v]], before[[v]], info = v)
  }
})

test_that("shut-down via the negative-assimilation exit also zeroes gas exchange", {
  # The second exit (assim_max_ < 0) does call E_from_Soil_to_Root_Collar, so
  # E_up_ and soil_consumption_ are evaluated at the zero-uptake collar
  # potential. The leaf-side pair is assigned nowhere on that path, though.
  ind <- TF24_Individual(TF24_Model())
  ind$set_state("height", 5)
  env <- TF24_Environment()
  env$set_soil_number_of_depths(5)

  invisible(tf24_aux_at(ind, env, 0.30))   # populate first
  env$extrinsic_drivers_set_constant("SWdown", 0)
  dark <- tf24_aux_at(ind, env, 0.30)

  expect_equal(dark[["transpiration"]], 0)
  expect_equal(dark[["stom_cond_CO2"]], 0)
  # Computed rather than assigned on this path, so it lands at numerical zero.
  expect_lt(abs(dark[["E_up_"]]), 1e-10)
})

test_that("no leaf output is NA once the plant has been evaluated", {
  # NA transpiration silently poisons any stand-level Et sum, which is the
  # headline water flux for site-level work.
  ind <- TF24_Individual(TF24_Model())
  ind$set_state("height", 5)
  env <- TF24_Environment()
  env$set_soil_number_of_depths(5)

  for (theta in c(0.30, 0.02, 0.30, 0.05)) {
    aux <- tf24_aux_at(ind, env, theta)
    for (v in c("transpiration", "E_up_", "stom_cond_CO2", "assimilation",
                "profit")) {
      expect_false(is.na(aux[[v]]),
                   info = sprintf("%s at theta = %g", v, theta))
    }
  }
})

test_that("a shut-down plant stops depleting the soil", {
  # The behavioural consequence, at patch scale: with no rain and soil dry
  # enough to shut the plants down, cumulative resource depletion must stop
  # advancing rather than continue at the last wet-step rate.
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.02, 5))
  env$extrinsic_drivers_set_constant("rainfall", 0)

  p <- testing_parameters("TF24")
  p$max_patch_lifetime <- 10
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  out <- tidy(run(p, env, record = "steps"))

  flux <- out$env$soil_moist_cumulative_flux
  expect_equal(max(flux$sum_resource_depletion), 0)
  # And the column holds its water: nothing drains, nothing is taken up.
  expect_equal(min(out$env$soil_moist$soil_moist), 0.02, tolerance = 1e-6)
})

test_that("a leaf in the canopy's shade is solved at the light it stands in", {
  # Darkness is an attainable state and the leaf is solved there. With the
  # energy balance off, at zero light it fixes no carbon, pays its dark
  # respiration and shuts, on every shading model; and as the light falls
  # towards zero the answer goes continuously to that one.
  at <- function(light, shading, theta = 0.30, energy_balance = FALSE) {
    m <- TF24_Model()
    m$control$shading_model <- shading
    if (energy_balance) m$pars$use_energy_balance <- 1
    ind <- TF24_Individual(m)
    ind$set_state("height", 5)
    env <- TF24_Environment()
    env$set_soil_number_of_depths(5)
    env$set_fixed_environment(light, 150)
    tf24_aux_at(ind, env, theta)
  }
  R_d_25 <- TF24_Model()$pars$R_d_25
  for (shading in c("mean-light", "crown-centre", "deep-crown")) {
    for (theta in c(0.30, 0.10)) {
      dark <- at(0, shading, theta)
      label <- sprintf("%s at theta = %g", shading, theta)
      expect_true(all(is.finite(dark)), label = label)
      expect_equal(dark[["assimilation"]], -R_d_25, tolerance = 1e-12,
                   label = label)
      # A trace, the solve's own tolerance on a shut leaf: 1.1e-12 with the
      # stem's head on, 3e-8 of the flow at a light of 1.
      expect_lt(dark[["transpiration"]], 1e-11)
      expect_lt(dark[["stom_cond_CO2"]], 1e-8)
    }
    # Below an openness of 1e-4 the light still reaches the leaf.
    dim <- vapply(c(1e-6, 1e-5, 1e-4), function(l)
      at(l, shading)[["assimilation"]], 0)
    expect_true(all(diff(dim) > 0), label = shading)
    expect_equal(at(1e-12, shading)[["assimilation"]], -R_d_25,
                 tolerance = 1e-9, label = shading)
  }

  # With the energy balance on, the leaf is cooler than the air in the dark, so
  # its respiration is below R_d_25, and in wet soil the stomata stay a little
  # open. Asserted here is that it answers: every output finite, carbon lost
  # rather than gained, and assimilation rising with the light from zero.
  for (theta in c(0.30, 0.10)) {
    lights <- c(0, 1e-6, 1e-5, 1e-4)
    eb <- lapply(lights, at, shading = "mean-light", theta = theta,
                 energy_balance = TRUE)
    label <- sprintf("energy balance at theta = %g", theta)
    expect_true(all(vapply(eb, function(x) all(is.finite(x)), TRUE)),
                label = label)
    A <- vapply(eb, `[[`, 0, "assimilation")
    expect_lt(A[[1]], 0)
    expect_true(all(diff(A) > 0), label = label)
  }
})
