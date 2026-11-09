# With the energy balance on, light heats the leaf through two bands: the
# absorbed PAR and the near-infrared, which TF24 seats as a fixed multiple of
# the intercepted PPFD. The leaf's derivative surface carries both, NIR through
# its own slot in the parameter pack, so a column that moves every cohort's
# light (k_I) is the pinned-step central difference. Filled without the NIR
# slot, the k_I column misses NIR's share of the heat load.

eb_stand_parameters <- function() {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- 6
  lma <- p$strategy_default$pars[["lma"]]
  p <- add_strategies(p, trait_matrix(c(0, 6, lma),
                                      c("TF24_floor_lambda_o", "hmat", "lma")),
                      hyperpar = TF24_hyperpar, birth_rate = 10)
  p$strategies[[1]]$pars$use_energy_balance <- 1
  full <- p$node_schedule_times[[1]]
  p$node_schedule_times <-
    list(full[round(seq(1, length(full), length.out = 12))])
  p
}

eb_stand_run <- function(p, record = TRUE) {
  env <- Environment("TF24")
  env$set_soil_water_state(rep(0.428 * 0.5, env$get_soil_number_of_depths()))
  env$extrinsic_drivers_set_constant("rainfall", 1)
  run_scm(p, env,
          Control(node_density_coordinate = "birth_date", ode_tol_rel = 1e-6,
                  ode_tol_abs = 1e-6),
          keep_states = record)
}

test_that("with the energy balance on, the light field's column is the pinned-step difference", {
  p <- eb_stand_parameters()
  scm <- eb_stand_run(p)
  off <- p
  off$strategies[[1]]$pars$use_energy_balance <- 0
  # The balance is live here: it moves fitness by orders of magnitude.
  expect_gt(stand_fitness(eb_stand_run(off, record = FALSE))[[1]],
            10 * stand_fitness(scm)[[1]])

  metrics <- c("1.net_reproduction_ratio", "1.offspring_production")
  names <- c("k_I", "vcmax_25")
  g <- stand_gradient(scm, metrics = metrics, traits = paste0("1.", names),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))

  pinned <- scm$parameters
  pinned$ode_times <- scm$ode_times
  pinned$ode_step_sizes <- scm$ode_step_sizes
  for (name in names) {
    x <- pinned$strategies[[1]]$pars[[name]]
    h <- 1e-5 * abs(x)
    up <- down <- pinned
    up$strategies[[1]]$pars[[name]] <- x + h
    down$strategies[[1]]$pars[[name]] <- x - h
    fd <- (stand_fitness(eb_stand_run(up, record = FALSE)) -
           stand_fitness(eb_stand_run(down, record = FALSE))) / (2 * h)
    for (metric in metrics) {
      expect_lt(abs(g$gradient[metric, paste0("1.", name)] - fd[[metric]]),
                1e-5 * abs(fd[[metric]]), label = paste(metric, name))
    }
  }
})
