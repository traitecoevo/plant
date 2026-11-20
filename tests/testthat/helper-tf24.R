# TF24 stands for the gradient and invasion tests.

# The TF24 environment with every soil layer at half its saturated water
# content (0.428) and a constant rain of 1.
tf24_half_wet_environment <- function() {
  env <- Environment("TF24")
  env$set_soil_water_state(rep(0.428 * 0.5, env$get_soil_number_of_depths()))
  env$extrinsic_drivers_set_constant("rainfall", 1)
  env
}

# A TF24 stand that reproduces, for the fitness referees: one species at the
# floor price's lambda_o = 0, with hmat lowered to 6 m so cohorts mature inside
# its six years, dense enough at a birth rate of 10 to hold an equilibrium, on
# twelve introductions. `lma_factor` scales the default lma, and
# `energy_balance` switches the leaf's energy balance on.
tf24_fitness_stand_parameters <- function(birth_rate = 10, lma_factor = 1,
                                          energy_balance = FALSE) {
  p <- testing_parameters("TF24")
  p$max_patch_lifetime <- 6
  lma <- p$strategy_default$pars[["lma"]]
  p <- add_strategies(p, trait_matrix(c(0, 6, lma * lma_factor),
                                      c("TF24_floor_lambda_o", "hmat", "lma")),
                      hyperpar = default_hyperpar("TF24"), birth_rate = birth_rate)
  if (energy_balance) {
    p$strategies[[1]]$pars$use_energy_balance <- 1
  }
  thin_node_schedule(p, 12)
}

# `p` run in the half-wet environment at an ODE tolerance of 1e-6, with its
# states recorded unless `record` is FALSE.
tf24_fitness_stand_run <- function(p, record = TRUE) {
  run(p, tf24_half_wet_environment(),
      Control(node_density_coordinate = "birth_date", ode_tol_rel = 1e-6,
              ode_tol_abs = 1e-6),
      record = if (record) "states" else character(0))$solver
}
