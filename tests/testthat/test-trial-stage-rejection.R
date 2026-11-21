# A Runge-Kutta step evaluates rates at trial stages the run will never keep.
# A stage can carry a newborn cohort to a density, or an individual to a size,
# the model has no meaning for. A guard reached there throws odelia's domain
# error, which rejects the step so the stepper retries it smaller; the same
# guard reached where the step cannot shrink, or outside a step, is an error
# that carries the guard's message.

# TF24 at its defaults with the seed mass alone lowered to 1.2 mg. On the first
# step after the introduction at t = 3, a stage takes the newborn's storage
# negative and its log density to hundreds, so its contribution to the
# competition integral is not finite; every accepted state stays finite.
tiny_seed_parameters <- function() {
  p <- testing_parameters("TF24")
  p$max_patch_lifetime <- 10
  add_strategies(p, trait_matrix(1.2e-6, "omega"), hyperpar = NULL)
}

test_that("a trial stage that overshoots a newborn rejects the step, and the run finishes", {
  skip_on_cran()
  p <- tiny_seed_parameters()
  r <- run(p)
  expect_true(r$complete)
  fitness <- r$net_reproduction_ratios[[1]]
  expect_true(is.finite(fitness))
  expect_gt(fitness, 0)
  # The same stand at a tolerance tight enough that no stage overshoots: the
  # rejected steps change which steps the run takes, nothing else.
  tight <- run(p, control = control(ode_tol_rel = 1e-7, ode_tol_abs = 1e-7))
  expect_equal(fitness, tight$net_reproduction_ratios[[1]], tolerance = 1e-2)
})

test_that("a guard reached outside a step is an error that names the state", {
  for (model in c("FF16", "TF24")) {
    env_name <- if (model == "FF16") "FF16_Env" else "TF24_Env"
    s <- get_list_of_strategy_types()[[model]]()
    s$birth_rate_y <- 1
    s$is_variable_birth_rate <- FALSE
    p <- Parameters(model, env_name)(strategies = list(s),
                                     fitness_weighting = "patch_age")
    names_ode <- Node(model, env_name)(s)$ode_names
    for (coord in c("height", "birth_date")) {
      info <- sprintf("%s on the %s coordinate", model, coord)
      patch <- Patch(model, env_name)(p, Environment(model),
                                      Control(node_density_coordinate = coord))
      patch$introduce_new_node(1, 0)
      patch$introduce_new_node(1, 1)
      y0 <- patch$ode_state

      # A density above the ceiling: the state check refuses it before any
      # rate is taken, on either coordinate.
      y <- y0
      y[[which(names_ode == "log_density")]] <- 60
      expect_error(patch$set_ode_state(y, 2),
                   "^Non-finite cohort density .* species 1 has a node with density=",
                   class = "odelia::util::DomainError", info = info)

      # A density inside the ceiling on an individual of infinite size: the
      # competition integral refuses the node's contribution, naming the node.
      y <- y0
      y[[which(names_ode == "height")]] <- Inf
      expect_error(patch$set_ode_state(y, 2),
                   "^Non-finite contribution to the competition integral: node 1 of 2 \\(height inf m",
                   class = "odelia::util::DomainError", info = info)
    }
  }
})
