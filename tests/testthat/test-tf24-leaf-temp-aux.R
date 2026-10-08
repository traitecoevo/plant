# The `Tleaf` auxiliary variable (#625).
#
# phylloptim's Leaf carries the leaf's own temperature at the operating point as
# an OUTPUT (`Tleaf_`), distinct from the `leaf_temp` driver it is handed. TF24
# computed it, used it to re-derive the whole Farquhar temperature block, and
# then discarded it -- so a canopy-level analysis could see the assimilation the
# temperature produced but not the temperature itself. These tests pin the aux
# slot's content, its units, and its behaviour under each shading model.
#
# Units: deg C.
#
# The distinction the slot exists for only becomes visible with
# pars$use_energy_balance non-zero. With it off (TF24's default) the leaf runs at
# the prescribed driver, and `Tleaf` is that driver exactly -- reported anyway,
# rather than NA, so the column can be plotted against anything.

tf24_temp_aux <- function(strategy = TF24_Strategy(), height = 5, ppfd = 1800,
                          theta = rep(0.25, 5), leaf_temp = NULL,
                          light = NULL) {
  ind <- TF24_Individual(strategy)
  ind$set_state("height", height)
  env <- TF24_Environment()
  env$extrinsic_drivers_set_constant("PPFD", ppfd)
  if (!is.null(leaf_temp)) {
    env$extrinsic_drivers_set_constant("leaf_temp", leaf_temp)
  }
  # `light` is canopy openness: a scalar for a uniform environment, or a
  # function of height for a vertical gradient (which is the only way the three
  # shading models can disagree -- see the deep-crown test below).
  if (is.function(light)) {
    hh <- seq(0, height, length.out = 101)
    ip <- Interpolator()
    ip$init(hh, vapply(hh, light, numeric(1)))
    env$light_availability$spline <- ip
  } else if (!is.null(light)) {
    env$set_fixed_environment(light, height_max = height)
  }
  env$set_soil_number_of_depths(length(theta))
  env$set_soil_water_state(theta)
  ind$compute_rates(env)
  stats::setNames(ind$internals$auxs, ind$aux_names)
}

pm_strategy <- function() {
  s <- TF24_Strategy()
  s$pars$use_energy_balance <- 1
  s
}

test_that("Tleaf aux is written, not left unset", {
  aux <- tf24_temp_aux()
  expect_true("Tleaf" %in% names(aux))
  expect_true(is.finite(aux[["Tleaf"]]))
})

test_that("with the energy balance off, Tleaf is the prescribed driver", {
  # Not a tautology worth skipping: it is the assertion that the slot tracks the
  # temperature the leaf actually ran at, and it is the only branch where the
  # expected value is known exactly.
  for (temp in c(15, 25, 35)) {
    aux <- tf24_temp_aux(leaf_temp = temp)
    expect_equal(aux[["Tleaf"]], temp, tolerance = 1e-12,
                 info = paste("leaf_temp =", temp))
  }
})

test_that("with the energy balance on, Tleaf departs from air temperature", {
  # The driver is reinterpreted as AIR temperature on this path, and a bright,
  # well-watered leaf sits above it. This is what the aux exists to report: off
  # the PM path the two are equal by construction, so a test that only ran with
  # the balance off could not distinguish an output from an echo of the input.
  air <- 30
  hot <- tf24_temp_aux(strategy = pm_strategy(), ppfd = 1800, leaf_temp = air)
  expect_true(is.finite(hot[["Tleaf"]]))
  expect_gt(hot[["Tleaf"]], air)

  # ...and a dark leaf is not warmed by radiation it does not absorb, so it sits
  # at or below air temperature. The pair brackets the driver, which is the
  # property that would fail if the slot were echoing `leaf_temp`.
  dark <- tf24_temp_aux(strategy = pm_strategy(), ppfd = 0, leaf_temp = air)
  expect_lte(dark[["Tleaf"]], hot[["Tleaf"]])
  expect_gt(hot[["Tleaf"]] - dark[["Tleaf"]], 1)
})

test_that("Tleaf rises with absorbed radiation under the energy balance", {
  by_light <- vapply(c(0, 200, 800, 1800),
                     function(q) tf24_temp_aux(strategy = pm_strategy(),
                                               ppfd = q)[["Tleaf"]],
                     numeric(1))
  expect_true(all(is.finite(by_light)))
  expect_true(all(diff(by_light) > 0))
})

test_that("Tleaf is finite when the hydraulic solve shuts down", {
  # Both shut-down exits bypass the transpiring branch, so they seat Tleaf at the
  # zero-transpiration (hottest) temperature explicitly. Dry top layers over wet
  # deep ones drives the shutdown; the same driver as the assimilation-aux test.
  dry <- c(rep(0.04, 8), rep(0.35, 7))
  aux <- tf24_temp_aux(theta = dry)
  expect_true(is.finite(aux[["Tleaf"]]))

  pm <- tf24_temp_aux(strategy = pm_strategy(), theta = dry, leaf_temp = 30)
  expect_true(is.finite(pm[["Tleaf"]]))
  # No latent cooling at zero transpiration, so the shut-down leaf is warmer than
  # air rather than pinned to it.
  expect_gt(pm[["Tleaf"]], 30)
})

test_that("Tleaf is reported under every shading model", {
  # Deep-crown must integrate Tleaf to a leaf-area-weighted crown mean alongside
  # the other leaf outputs. Left out of that loop it would report whichever
  # quadrature node the loop ended on -- finite, plausible, and wrong.
  for (model in c("mean-light", "crown-centre", "deep-crown")) {
    s <- pm_strategy()
    s$control$shading_model <- model
    aux <- tf24_temp_aux(strategy = s, leaf_temp = 30)
    expect_true(is.finite(aux[["Tleaf"]]), info = model)
    expect_gt(aux[["Tleaf"]], 30)
  }
})

test_that("deep-crown Tleaf is a crown mean over the light gradient", {
  # Under a UNIFORM light environment every quadrature node sees the same light,
  # so all three models agree exactly and this test could not tell an integral
  # from a single node. The gradient is what makes the assertion mean something:
  # a stand-alone Individual sits in no patch, so it is imposed directly on the
  # environment's light spline.
  h <- 15
  light <- function(z) 0.08 + 0.9 * (z / h)^2   # dark base, bright top
  make <- function(model, ...) {
    s <- pm_strategy()
    s$control$shading_model <- model
    tf24_temp_aux(strategy = s, height = h, leaf_temp = 30, ...)[["Tleaf"]]
  }

  deep <- make("deep-crown", light = light)
  # The coolest and warmest leaves in the crown: the same plant under uniform
  # light equal to the crown base's and the crown top's. A leaf-area-weighted
  # mean of the profile must sit strictly inside that bracket -- which is exactly
  # what leaving Tleaf at one quadrature node would not guarantee.
  base <- make("deep-crown", light = light(0))
  top <- make("deep-crown", light = light(h))
  expect_true(all(is.finite(c(deep, base, top))))
  expect_gt(top, base)
  expect_gt(deep, base)
  expect_lt(deep, top)

  # ...and the gradient makes the three shading models genuinely disagree, so a
  # deep-crown value equal to either single-evaluation model to the last digit
  # would mean the integral had collapsed onto one node.
  expect_false(isTRUE(all.equal(deep, make("crown-centre", light = light),
                                tolerance = 1e-10)))
  expect_false(isTRUE(all.equal(deep, make("mean-light", light = light),
                                tolerance = 1e-10)))
})

test_that("under uniform light every shading model gives the same Tleaf", {
  # The counterpart of the test above, and the property that pins the weighting:
  # q integrates to one over the crown, so the deep-crown integral of a constant
  # is that constant. An un-normalised sum would not be.
  make <- function(model) {
    s <- pm_strategy()
    s$control$shading_model <- model
    tf24_temp_aux(strategy = s, height = 15, leaf_temp = 30,
                  light = 0.4)[["Tleaf"]]
  }
  ref <- make("deep-crown")
  expect_equal(make("crown-centre"), ref, tolerance = 1e-10)
  expect_equal(make("mean-light"), ref, tolerance = 1e-10)
})
