# Fitness as a differentiable stand metric on TF24: each species' net
# reproduction ratio and offspring production, refereed against a central
# difference of runs pinned to the recorded steps, as FF16's are
# (test-ff16-fitness-gradient.R). Each row is refereed on its own scale: the
# fitness sits orders of magnitude away from the census.
#
# The stand reproduces: hmat is lowered to 6 m so cohorts mature inside its six
# years, and it is dense enough at a birth rate of 10 to hold an equilibrium.

fitness_stand_parameters <- function(birth_rate, lma_factor = 1) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- 6
  lma <- p$strategy_default$pars[["lma"]]
  p <- add_strategies(p, trait_matrix(c(0, 6, lma * lma_factor),
                                      c("TF24_floor_lambda_o", "hmat", "lma")),
                      hyperpar = TF24_hyperpar, birth_rate = birth_rate)
  full <- p$node_schedule_times[[1]]
  p$node_schedule_times <-
    list(full[round(seq(1, length(full), length.out = 12))])
  p
}

fitness_stand_environment <- function() {
  env <- Environment("TF24")
  env$set_soil_water_state(rep(0.428 * 0.5, env$get_soil_number_of_depths()))
  env$extrinsic_drivers_set_constant("rainfall", 1)
  env
}

fitness_stand_control <- function() {
  Control(node_density_in_birth_date = TRUE, ode_tol_rel = 1e-6,
          ode_tol_abs = 1e-6)
}

fitness_stand_run <- function(p, record = TRUE) {
  run_scm(p, fitness_stand_environment(), fitness_stand_control(),
          record_trajectory = record)
}

test_that("TF24 reports its fitness under the names the gradient uses", {
  scm <- fitness_stand_run(fitness_stand_parameters(10))
  f <- stand_fitness(scm)
  expect_equal(names(f), c("1.net_reproduction_ratio", "1.offspring_production"))
  expect_identical(f[["1.net_reproduction_ratio"]], scm$net_reproduction_ratios)
  expect_identical(f[["1.offspring_production"]], scm$offspring_production)
  expect_false("S_D" %in% names(census_undifferentiable_tf24()))
  g <- stand_gradient(scm, metrics = names(f), hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  # The birth rate is the last column of the species' block.
  expect_equal(utils::tail(colnames(g$gradient), 1), "1.birth_rate")
  expect_true("1.S_D" %in% colnames(g$gradient))
})

test_that("TF24's fitness gradient is the pinned-step central difference", {
  scm <- fitness_stand_run(fitness_stand_parameters(10))
  metrics <- c("1.net_reproduction_ratio", "1.offspring_production")
  # Every route into the fecundity state: the reproduction schedule (hmat,
  # a_f1, a_f2), the per-seed cost (omega, a_f3), survival during dispersal,
  # growth and the light field (lma, k_I), mortality (d_I, a_dG1), the leaf
  # (vcmax_25, stem_P50), establishment (a_d0) and the birth rate.
  names <- c("hmat", "a_f1", "a_f2", "omega", "a_f3", "S_D", "lma", "k_I",
             "d_I", "a_dG1", "vcmax_25", "stem_P50", "a_d0")
  g <- stand_gradient(scm, metrics = metrics,
                      traits = paste0("1.", c(names, "birth_rate")),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  expect_gt(g$value[["1.net_reproduction_ratio"]], 1)

  pinned <- scm$parameters
  pinned$ode_times <- scm$ode_times
  pinned$ode_step_sizes <- scm$ode_step_sizes
  # lma and omega set the seed height, which is solved to a root tolerance that
  # a difference reads: omega's is off by 2.5e-2 at a step of 1e-5 and 8e-4 at
  # 1e-3. Every other column is differenced at 1e-4 and held to 1e-4.
  birth_size <- c("lma", "omega")
  difference <- function(set, x, rel) {
    h <- rel * abs(x)
    up <- stand_fitness(fitness_stand_run(set(x + h), record = FALSE))
    down <- stand_fitness(fitness_stand_run(set(x - h), record = FALSE))
    (up - down) / (2 * h)
  }
  with_parameter <- function(name) function(v) {
    q <- pinned
    q$strategies[[1]]$pars[[name]] <- v
    q
  }
  with_birth_rate <- function(v) {
    q <- pinned
    q$strategies[[1]]$birth_rate_y <- v
    q
  }
  for (name in c(names, "birth_rate")) {
    rel <- if (name %in% birth_size) 1e-3 else 1e-4
    fd <- if (name == "birth_rate") {
      difference(with_birth_rate, pinned$strategies[[1]]$birth_rate_y, rel)
    } else {
      difference(with_parameter(name), pinned$strategies[[1]]$pars[[name]], rel)
    }
    tolerance <- if (name %in% birth_size) 5e-3 else 1e-4
    for (metric in metrics) {
      expect_lt(abs(g$gradient[metric, paste0("1.", name)] - fd[[metric]]),
                tolerance * abs(fd[[metric]]), label = paste(metric, name))
    }
  }
  # S_D multiplies the whole integral, so its column is the value over S_D.
  S_D <- pinned$strategies[[1]]$pars$S_D
  expect_equal(g$gradient["1.net_reproduction_ratio", "1.S_D"],
               g$value[["1.net_reproduction_ratio"]] / S_D, tolerance = 1e-12)
})

test_that("Newton on d ln J / d b finds the equilibrium birth rate, and an invader is differentiated there", {
  # The resident at its demographic equilibrium makes as many seeds as it is
  # given: J(b*) = b*. Newton in ln b on ln J - ln b, with the slope from the
  # birth-rate column of offspring production J.
  b <- 10
  for (k in 1:8) {
    resident <- fitness_stand_run(fitness_stand_parameters(b))
    g <- stand_gradient(resident, metrics = "1.offspring_production",
                        traits = "1.birth_rate", hyperpar = NULL)
    J <- g$value[["1.offspring_production"]]
    excess <- log(J) - log(b)
    if (abs(excess) < 1e-8) {
      break
    }
    b <- exp(log(b) - excess / (b * g$gradient[[1]] / J - 1))
  }
  # Linear, not quadratic: the adaptive runs' slope in b differs from the
  # swept slope of their own steps by about 2.5%. Seven runs here.
  expect_lt(abs(excess), 1e-8)
  expect_equal(resident$net_reproduction_ratios, 1, tolerance = 1e-8)

  # An invader at the resident's traits has the resident's ratio, one.
  resident$record_trajectory <- TRUE
  resident$run_mutant(fitness_stand_parameters(1))
  expect_equal(resident$net_reproduction_ratios, 1, tolerance = 1e-8)

  # An invader at lma x0.9, differentiated with the canopy held at the
  # resident's equilibrium: the sweep against a central difference of
  # invasions of the same recording.
  invader <- fitness_stand_parameters(1, 0.9)
  resident$run_mutant(invader)
  gi <- stand_gradient(resident, metrics = "1.net_reproduction_ratio",
                       traits = c("1.lma", "1.hmat"), hyperpar = NULL)
  expect_false(any(stand_gradient_refused(gi)))
  expect_gt(gi$value[["1.net_reproduction_ratio"]], 1)
  resident$record_trajectory <- FALSE
  for (name in c("lma", "hmat")) {
    x <- invader$strategies[[1]]$pars[[name]]
    h <- 1e-4 * x
    at <- function(v) {
      q <- invader
      q$strategies[[1]]$pars[[name]] <- v
      resident$run_mutant(q)
      resident$net_reproduction_ratios
    }
    fd <- (at(x + h) - at(x - h)) / (2 * h)
    expect_equal(gi$gradient["1.net_reproduction_ratio", paste0("1.", name)],
                 fd, tolerance = if (name == "lma") 1e-3 else 1e-4, label = name)
  }
})
