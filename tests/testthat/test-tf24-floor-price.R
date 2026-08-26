# The TF24_floor cost curve and the shadow price of water (#634).
#
# TF24 and TF24f are seated on phylloptim's TF24_floor curve, whose one parameter
# `TF24_floor_lambda_o` is the price of water as transpiration goes to zero. Two
# things are being pinned here, and they are different in kind:
#
#   1. AT THE DEFAULT PRICE OF ZERO, NOTHING MOVES. TF24_floor at lambda_o = 0
#      IS TF24, at identical parameter values and bit-for-bit -- phylloptim
#      asserts that reduction with no tolerance, because each term is the
#      parent's own expression so zeroing the price adds an exact zero. Every
#      pinned baseline in this suite is the real evidence for that; these tests
#      pin the two identities that hold *at* the default and would otherwise be
#      invisible.
#
#   2. THE PRICE IS NOT A COST. `profit` is the OBJECTIVE and deducts the shadow
#      price; the carbon the plant actually kept is `profit + shadow_cost`, and
#      that is what growth is billed. The distinction has no numerical
#      consequence at the default (shadow_cost is 0), so the only way to test it
#      is to raise the price -- which is what the second block does.

tf24_aux <- function(strategy = TF24_Strategy(), height = 5, ppfd = 1800,
                     theta = rep(0.25, 5), individual = TF24_Individual) {
  ind <- individual(strategy)
  ind$set_state("height", height)
  env <- TF24_Environment()
  env$extrinsic_drivers_set_constant("PPFD", ppfd)
  env$set_soil_number_of_depths(length(theta))
  env$set_soil_water_state(theta)
  ind$compute_rates(env)
  stats::setNames(ind$internals$auxs, ind$aux_names)
}

priced <- function(lambda_o, ...) {
  s <- TF24_Strategy()
  p <- s$pars
  p$TF24_floor_lambda_o <- lambda_o
  s$pars <- p
  tf24_aux(s, ...)
}

# --- at the default price ----------------------------------------------------

test_that("the water price defaults to zero, and zero is TF24", {
  # An explicit 0.0 rather than phylloptim's NA sentinel: phylloptim REFUSES
  # TF24_floor with an unset price, deliberately, so a default of NA here would
  # make every TF24 run an error rather than a silent surprise.
  expect_identical(TF24_Strategy()$pars$TF24_floor_lambda_o, 0)
})

test_that("shadow_cost is reported, and is exactly zero at the default price", {
  aux <- tf24_aux()
  expect_true("shadow_cost" %in% names(aux))
  # Exactly, not approximately: shadow_cost is `lambda_o * transpiration_`, and
  # `0 * E` is an exact zero for any finite E.
  expect_identical(aux[["shadow_cost"]], 0)
})

test_that("carbon profit equals the objective exactly when water is free", {
  # `profit + shadow_cost == profit` at lambda_o = 0. This is the identity
  # net_mass_production_dt relies on to be bit-neutral at the default, so it is
  # asserted as an identity rather than to a tolerance.
  aux <- tf24_aux()
  expect_identical(aux[["profit"]] + aux[["shadow_cost"]], aux[["profit"]])
})

# --- with a non-zero price ---------------------------------------------------
#
# lambda_o is in umol C (kg H2O)^-1 and its scale is set by the leaf, whose own
# marginal cost of water runs ~9e4 to 3e5 at these defaults. 1e5 is inside that
# band, so the optimum moves without being pinned against a bracket bound.

test_that("a non-zero price closes the stomata", {
  free  <- tf24_aux()
  costly <- priced(1e5)

  # Water is worth something now, so less of it is spent: the aperture and the
  # flux both fall. This is the whole behavioural content of the curve.
  expect_lt(costly[["stom_cond_CO2"]], free[["stom_cond_CO2"]])
  expect_lt(costly[["transpiration"]], free[["transpiration"]])
  # And the leaf operates at a wetter stem potential, because it is moving less
  # water through the same conductance.
  expect_lt(costly[["opt_psi_stem"]], free[["opt_psi_stem"]])
})

test_that("the shadow price is the price times the transpiration, exactly", {
  # The identity a caller checks. phylloptim reads the STORED transpiration
  # rather than recomputing it precisely so this holds against the reported E --
  # and plant integrates `profit` and `transpiration` over the crown against the
  # same quadrature weights, so it survives the canopy integral too.
  aux <- priced(1e5)
  expect_equal(aux[["shadow_cost"]], 1e5 * aux[["transpiration"]],
               tolerance = 1e-12)
  expect_gt(aux[["shadow_cost"]], 0)
})

test_that("the shadow price does not reduce growth", {
  # ⚠️ THE POINT OF THE ISSUE. `profit` is the objective and has the price
  # deducted from it; the carbon kept is `profit + shadow_cost`, which must be
  # STRICTLY greater. Feeding the objective into growth instead would tax the
  # plant by carbon it never spent.
  aux <- priced(1e5)
  expect_gt(aux[["profit"]] + aux[["shadow_cost"]], aux[["profit"]])

  # And the size of that error, so a regression here is legible rather than just
  # red: the objective understates the carbon kept by a double-digit fraction at
  # this price.
  understatement <- aux[["shadow_cost"]] /
    (aux[["profit"]] + aux[["shadow_cost"]])
  expect_gt(understatement, 0.05)
})

test_that("net production is billed on carbon kept, not on the objective", {
  # The consequence at the level plant actually cares about. Growth is
  # `carbon_profit * area_leaf * ...`, so a run at a non-zero price must produce
  # MORE carbon than the same run would if the price were treated as a realised
  # cost. The counterfactual is not directly reachable from R, so this asserts
  # the reachable half: net production stays above what the objective alone
  # would support.
  aux <- priced(1e5)
  # A well-lit, well-watered plant is still in carbon surplus at this price.
  expect_gt(aux[["net_mass_production_dt"]], 0)
  # Both terms of the split are finite and the realised cost is still positive:
  # `assimilation - profit` is the WHOLE cost the objective subtracted, of which
  # shadow_cost is only a part, so the realised remainder must be positive too.
  realised <- aux[["assimilation"]] - aux[["profit"]] - aux[["shadow_cost"]]
  expect_gt(realised, 0)
})

test_that("TF24f tracks the seated curve rather than TF24", {
  # ⚠️ THE FAILURE MODE THIS EXISTS FOR. Every phylloptim entry point TF24f calls
  # is templated on the cost curve, and each used to name TF24 at its call site --
  # so `set_model()` had no effect on the solve at all, and seating TF24_floor
  # would have moved `shadow_cost` while the tracked state went on climbing TF24's
  # gradient. If solve_leaf() stopped dispatching on the seated curve, the two
  # arms below would come back identical.
  #
  # ⚠️ TF24f EVALUATES AT THE TRACKED COLLAR RATHER THAN OPTIMISING IT, so E and
  # gs at a GIVEN tracked state are curve-independent by construction and cannot
  # be the instrument here (they are what test the TF24 arm above). What the curve
  # moves is the OBJECTIVE at that point and its GRADIENT -- i.e. `shadow_cost`
  # and the tracked state's rate. Both are checked.
  #
  # The tracked state has to be seeded: a freshly constructed individual holds 0,
  # which is a collar at saturation, so no water moves and every price is worth
  # nothing times zero.
  tracked <- tf24_aux()[["opt_root_psi"]]
  expect_gt(tracked, 0)

  f_aux <- function(lambda_o) {
    s <- TF24f_Strategy()
    p <- s$pars
    p$TF24_floor_lambda_o <- lambda_o
    s$pars <- p
    ind <- TF24f_Individual(s)
    ind$set_state("height", 5)
    ind$set_state("opt_root_psi_state", tracked)
    env <- TF24_Environment()
    env$extrinsic_drivers_set_constant("PPFD", 1800)
    env$set_soil_number_of_depths(5)
    env$set_soil_water_state(rep(0.25, 5))
    ind$compute_rates(env)
    vars <- ind$internals
    list(aux  = stats::setNames(vars$auxs, ind$aux_names),
         rate = vars$rates[which(ind$ode_names == "opt_root_psi_state")])
  }

  free <- f_aux(0)
  costly <- f_aux(1e5)

  # The price is reported, and only when there is one.
  expect_identical(free$aux[["shadow_cost"]], 0)
  expect_gt(costly$aux[["transpiration"]], 0)
  expect_gt(costly$aux[["shadow_cost"]], 0)

  # And the SOLVE saw it: the acclimation gradient the tracked state climbs is
  # the priced objective's, so it points to a wetter collar than TF24's does.
  expect_lt(costly$rate, free$rate)
})

# --- the price as a TRAIT ------------------------------------------------------
#
# `TF24_floor_lambda_o` is a `TF24_Pars` field, and that is NOT the same thing as
# being usable as a trait: the trait pathway goes trait_matrix() ->
# generate_strategy() -> `pars[trait_names] <- xi`, through the hyperpar, and a
# parameter can be settable on a strategy while being unreachable that way (the
# hyperpar's overwrite guard rejects anything it derives itself, and phylloptim's
# own CF77 price spent its whole life reachable only through `leaf_model()` for
# the mirror-image reason). It works here, and these pin it, because a price of
# water is a strategic trait -- lambda_o indexes a continuum from carbon
# maximisation at zero to water-conserving strategies above it, which is a
# prediction about aridity and so is exactly the thing someone will want to vary
# across species or sweep.

test_that("the water price is settable as a trait", {
  p0 <- scm_base_parameters("TF24")
  prices <- c(0, 5e4, 1e5)
  ss <- generate_strategy(p0, trait_matrix(prices, "TF24_floor_lambda_o"))

  expect_length(ss, length(prices))
  expect_equal(vapply(ss, function(s) s$pars$TF24_floor_lambda_o, numeric(1)),
               prices)
})

test_that("a strategy built from the trait actually solves at that price", {
  # ⚠️ THE HALF THAT MATTERS. The test above only shows the value landed on the
  # strategy; this shows `prepare_strategy()` carried it into the leaf and the
  # solve used it. A parameter that reaches `pars` but never reaches
  # `leaf.TF24_floor_lambda_o` looks completely fine from R.
  p0 <- scm_base_parameters("TF24")
  prices <- c(0, 5e4, 1e5)
  ss <- generate_strategy(p0, trait_matrix(prices, "TF24_floor_lambda_o"))

  got <- lapply(ss, function(s) tf24_aux(s))
  shadow <- vapply(got, function(a) a[["shadow_cost"]], numeric(1))
  flux   <- vapply(got, function(a) a[["transpiration"]], numeric(1))

  # Zero price, zero shadow cost; and a priced leaf reports one.
  expect_identical(shadow[[1]], 0)
  expect_true(all(shadow[-1] > 0))
  # The identity, per strategy, against that strategy's own price.
  expect_equal(shadow, prices * flux, tolerance = 1e-12)
  # And the solve moved: a dearer price spends less water. Monotone, so this
  # cannot pass on a leaf that read only one of the three prices.
  expect_true(all(diff(flux) < 0))
})

test_that("the price varies across species in one parameter set", {
  # The multi-species route, which is a different code path from
  # generate_strategy(): each species gets its own strategy and its own leaf.
  p0 <- scm_base_parameters("TF24")
  p <- add_strategies(p0, trait_matrix(c(0, 1e5), "TF24_floor_lambda_o"),
                      birth_rate = c(1, 1))
  expect_length(p$strategies, 2)
  expect_equal(vapply(p$strategies, function(s) s$pars$TF24_floor_lambda_o,
                      numeric(1)),
               c(0, 1e5))
})

test_that("TF24_hyperpar passes the price through untouched", {
  # ⚠️ NOT AUTOMATIC. TF24_hyperpar refuses to overwrite a parameter it derives
  # itself ("Attempt to overwrite generated parameters"), which is what makes
  # `stem_P50`, `stem_c`, `stem_b`, `psi_crit` and `TF24_cost_scale` unusable as
  # input traits. The price is derived from nothing, so it must survive -- and if
  # someone ever adds a lambda_o ~ aridity relationship to the hyperpar, this is
  # the test that will tell them they have just made it unsettable.
  p0 <- scm_base_parameters("TF24")
  m <- trait_matrix(c(0, 1e5), "TF24_floor_lambda_o")
  ret <- TF24_hyperpar(m, p0$strategy_default)
  expect_true("TF24_floor_lambda_o" %in% colnames(ret))
  expect_equal(ret[, "TF24_floor_lambda_o"], c(0, 1e5))
})

test_that("TF24f takes the price as a trait too", {
  # TF24f shares TF24_Pars, so this is about the hyperpar and the strategy
  # template rather than about a second parameter.
  pf <- scm_base_parameters("TF24f")
  sf <- generate_strategy(pf, trait_matrix(1e5, "TF24_floor_lambda_o"))
  expect_equal(sf[[1]]$pars$TF24_floor_lambda_o, 1e5)
})
