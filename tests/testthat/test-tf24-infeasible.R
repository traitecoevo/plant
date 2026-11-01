# An infeasible leaf solve reaches odelia as a domain error (#608). phylloptim's
# infeasible_error is a std::runtime_error and not an odelia DomainError, so
# untranslated it ends a run having taken zero steps, where translated the
# stepper shrinks the step and retries. The translation says so in its message,
# which is what these check, for both solve_leaf bodies.

infeasible_rates <- function(strategy, individual) {
  # A dry, dark or hot leaf is a value rather than an infeasibility since
  # phylloptim 0.19, so what is left is a solver that cannot converge: a root-find
  # held to one iteration, which phylloptim reports as infeasible.
  ctrl <- strategy$control
  ctrl$ci_niter <- 1
  strategy$control <- ctrl
  ind <- individual(strategy)
  ind$set_state("height", 5)
  env <- TF24_Environment()
  env$extrinsic_drivers_set_constant("SWdown", 1800 / (0.5 * 4.57))  # PPFD above the canopy, as SWdown
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.25, 5))
  ind$compute_rates(env)
}

test_that("TF24 translates an infeasible leaf solve into a domain error", {
  expect_error(infeasible_rates(TF24_Strategy(), TF24_Individual),
               "^leaf solve infeasible: \\[phylloptim:infeasible:")
})

test_that("TF24f translates an infeasible leaf solve into a domain error", {
  expect_error(infeasible_rates(TF24f_Strategy(), TF24f_Individual),
               "^leaf solve infeasible: \\[phylloptim:infeasible:")
})
