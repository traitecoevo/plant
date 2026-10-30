# Fitness as a differentiable stand metric on FF16: each species' net
# reproduction ratio and offspring production, with the census's shape,
# refereed the way the census is -- a central difference of runs pinned to the
# recorded steps at a tight seed-height root tolerance (see
# test-ff16-stand-gradient.R for why both).
#
# The stand is built to reproduce: hmat is lowered to 3 m so cohorts mature
# inside the run, which is what puts the reproduction parameters on the path.

fitness_control <- function(...) {
  Control(node_density_coordinate = "birth_date", offspring_production_tol = 1e-12,
          ...)
}

fitness_parameters <- function(lifetime = 8, lma = 0.0825, hmat = 3,
                               birth_rate = 2) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  m <- trait_matrix(c(lma, rep(hmat, length(lma))), c("lma", "hmat"))
  add_strategies(p, m, birth_rate = as.list(rep(birth_rate, length(lma))))
}

fitness_stand <- function(..., ctrl = fitness_control()) {
  run_scm(fitness_parameters(...), NULL, ctrl, keep_states = TRUE)
}

pinned <- function(scm) {
  p <- scm$parameters
  p$ode_times <- scm$ode_times
  p$ode_step_sizes <- scm$ode_step_sizes
  p
}

# d(metric)/d(parameter) of species `species`, the perturbed runs stepping
# where the recorded one stepped.
central_difference <- function(p, species, name, metric, rel = 1e-4) {
  x <- p$strategies[[species]]$pars[[name]]
  h <- rel * abs(x)
  at <- function(v) {
    q <- p
    q$strategies[[species]]$pars[[name]] <- v
    stand_fitness(run_scm(q, NULL, fitness_control()))[[metric]]
  }
  (at(x + h) - at(x - h)) / (2 * h)
}

test_that("FF16 reports its fitness under the names the gradient uses", {
  scm <- fitness_stand(lma = c(0.0825, 0.2))
  f <- stand_fitness(scm)
  expect_equal(names(f), c("1.net_reproduction_ratio", "1.offspring_production",
                           "2.net_reproduction_ratio", "2.offspring_production"))
  # The same numbers the run has always reported, bit for bit.
  expect_identical(unname(f[c(1, 3)]), scm$net_reproduction_ratios)
  expect_identical(unname(f[c(2, 4)]), scm$offspring_production)
  # The census keeps its own three, and the gradient defaults to all seven.
  expect_equal(names(stand_census(scm)),
               c("leaf_area", "mass_above_ground", "area_stem"))
  g <- stand_gradient(scm, traits = "1.lma", hyperpar = NULL)
  expect_equal(rownames(g$gradient), c(names(stand_census(scm)), names(f)))
  expect_false(any(stand_gradient_refused(g)))
  # S_D reaches the fitness directly, so it is a column now.
  expect_true("1.S_D" %in% colnames(stand_gradient(scm, metrics = "1.offspring_production",
                                                   hyperpar = NULL)$gradient))
  expect_error(stand_gradient(scm, metrics = "3.offspring_production",
                              hyperpar = NULL),
               "Unknown metric")
})

test_that("the fitness gradient is the pinned-step central difference", {
  scm <- fitness_stand()
  expect_gt(stand_fitness(scm)[["1.net_reproduction_ratio"]], 0)
  g <- stand_gradient(scm, metrics = c("1.net_reproduction_ratio",
                                       "1.offspring_production"),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  p <- pinned(scm)
  # Every route into the fecundity state: the reproduction schedule (hmat,
  # a_f1, a_f2), the per-seed cost (omega, a_f3), survival during dispersal
  # (S_D, a direct term and nothing else), growth and the light field (lma,
  # k_I), mortality (a_dG1, d_I) and establishment (a_d0).
  for (name in c("hmat", "a_f1", "a_f2", "omega", "a_f3", "S_D", "lma", "k_I",
                 "a_dG1", "d_I", "a_d0")) {
    for (metric in c("1.net_reproduction_ratio", "1.offspring_production")) {
      fd <- central_difference(p, 1, name, metric)
      # Relative to the difference, not absolute: a_dG1's column is ~2.6e-6, so
      # an absolute tolerance of 1e-5 would pass zero or its negative.
      expect_lt(abs(g$gradient[metric, paste0("1.", name)] - fd), 1e-5 * abs(fd),
                label = paste(metric, name))
    }
  }
  # With a constant birth rate the offspring production is the net reproduction
  # ratio times that rate, column for column.
  # Every column but the birth rate's own, which carries the ratio itself too
  # (test-ff16-birth-rate-gradient.R has that one).
  b <- scm$parameters$strategies[[1]]$birth_rate_y
  traits <- setdiff(colnames(g$gradient), "1.birth_rate")
  expect_equal(g$gradient["1.offspring_production", traits],
               b * g$gradient["1.net_reproduction_ratio", traits],
               tolerance = 1e-12)
  # S_D multiplies the whole integral, so its column is the value over S_D.
  S_D <- scm$parameters$strategies[[1]]$pars$S_D
  expect_equal(g$gradient["1.net_reproduction_ratio", "1.S_D"],
               g$value[["1.net_reproduction_ratio"]] / S_D, tolerance = 1e-12)
})

test_that("one species' fitness responds to the other's traits through the canopy", {
  scm <- fitness_stand(lma = c(0.0825, 0.2))
  g <- stand_gradient(scm, metrics = c("1.net_reproduction_ratio",
                                       "2.net_reproduction_ratio"),
                      traits = c("1.lma", "2.lma", "1.k_I", "2.k_I"),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  p <- pinned(scm)
  # A finer step than the one-species test: one species' fitness is strongly
  # curved in the other's lma, and at a relative step of 1e-4 the difference's
  # own truncation error is 2e-5 (it falls as the square of the step: 3e-7 at
  # 1e-5), which is not the sweep's.
  for (i in 1:2) for (j in 1:2) for (name in c("lma", "k_I")) {
    metric <- paste0(i, ".net_reproduction_ratio")
    fd <- central_difference(p, j, name, metric, rel = 1e-5)
    expect_equal(g$gradient[metric, paste0(j, ".", name)], fd,
                 tolerance = 1e-5, label = paste(metric, j, name))
  }
  # The cross terms are real: species 2's shading reaches species 1's fitness.
  expect_true(abs(g$gradient["1.net_reproduction_ratio", "2.k_I"]) > 0)
})

test_that("the fitness gradient over a window is the pinned-step central difference", {
  # A window takes the weighted mean over its own cohorts. The weights are the
  # regime's, read at the birth dates, so they carry no derivative and the
  # normalisation divides the gradient as it divides the value.
  p <- fitness_window(fitness_parameters(), 2, 5)
  scm <- run_scm(p, NULL, fitness_control(), keep_states = TRUE)
  g <- stand_gradient(scm, metrics = c("1.net_reproduction_ratio",
                                       "1.offspring_production"),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  expect_identical(g$value[["1.net_reproduction_ratio"]],
                   scm$net_reproduction_ratios)
  q <- pinned(scm)
  for (name in c("hmat", "lma", "S_D")) {
    for (metric in c("1.net_reproduction_ratio", "1.offspring_production")) {
      fd <- central_difference(q, 1, name, metric)
      expect_equal(g$gradient[metric, paste0("1.", name)], fd,
                   tolerance = 1e-5, label = paste(metric, name))
    }
  }
  traits <- setdiff(colnames(g$gradient), "1.birth_rate")
  expect_equal(g$gradient["1.offspring_production", traits],
               2 * g$gradient["1.net_reproduction_ratio", traits],
               tolerance = 1e-12)
})
