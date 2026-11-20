# With the energy balance on, light heats the leaf through two bands: the
# absorbed PAR and the near-infrared, which TF24 seats as a fixed multiple of
# the intercepted PPFD. The leaf's derivative surface carries both, NIR through
# its own slot in the parameter pack, so a column that moves every cohort's
# light (k_I) is the pinned-step central difference. Filled without the NIR
# slot, the k_I column misses NIR's share of the heat load.

test_that("with the energy balance on, the light field's column is the pinned-step difference", {
  p <- tf24_fitness_stand_parameters(energy_balance = TRUE)
  solver <- tf24_fitness_stand_run(p)
  off <- p
  off$strategies[[1]]$pars$use_energy_balance <- 0
  # The balance is live here: it moves fitness by orders of magnitude.
  expect_gt(stand_fitness(tf24_fitness_stand_run(off, record = FALSE))[[1]],
            10 * stand_fitness(solver)[[1]])

  metrics <- c("1.net_reproduction_ratio", "1.offspring_production")
  names <- c("k_I", "vcmax_25")
  g <- stand_gradient(solver, metrics = metrics, traits = paste0("1.", names),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))

  pinned <- pin_steps(solver)
  fitness_of <- function(q) stand_fitness(tf24_fitness_stand_run(q, record = FALSE))
  for (name in names) {
    fd <- central_difference(pinned, 1, name, NULL, fitness_of, rel = 1e-5)
    for (metric in metrics) {
      expect_lt(abs(g$gradient[metric, paste0("1.", name)] - fd[[metric]]),
                1e-5 * abs(fd[[metric]]), label = paste(metric, name))
    }
  }
})
