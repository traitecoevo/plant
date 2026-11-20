# The birth-rate column on FF16: d(metric)/d(birth rate of
# species j), the Jacobian of the offspring map an equilibrium solver iterates
# on. A constant birth rate reaches the model at the inflow condition (the
# density every cohort is born at) and as the offspring integral's scalars, and
# is carried at the working scalar so both routes are on the tape. A variable
# birth rate is a driver and has no column.

# Two species, at birth rates `b`, on the reproducing stand.
birth_stand <- function(b) {
  ff16_stand_run(ff16_stand_parameters(c(0.0825, 0.2), birth_rate = b))
}

test_that("the birth rate is a column, named per species, for a constant rate", {
  solver <- birth_stand(c(2, 1))
  g <- stand_gradient(solver, metrics = "1.net_reproduction_ratio", hyperpar = NULL)
  expect_true(all(c("1.birth_rate", "2.birth_rate") %in% colnames(g$gradient)))
  expect_false(any(stand_gradient_refused(g)))
  # The last column of each species' block, after the parameter table's.
  cols <- colnames(g$gradient)
  expect_equal(which(cols == "1.birth_rate"), which(cols == "2.lma") - 1L)
})

test_that("d(fitness)/d(birth rate) is the pinned-step central difference", {
  b <- c(2, 1)
  solver <- birth_stand(b)
  metrics <- c("1.net_reproduction_ratio", "2.net_reproduction_ratio",
               "1.offspring_production", "2.offspring_production")
  g <- stand_gradient(solver, metrics = metrics,
                      traits = c("1.birth_rate", "2.birth_rate"),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  p <- pin_steps(solver)
  fitness_of <- function(q) stand_fitness(run(q, NULL, ff16_referee_control())$solver)
  for (j in 1:2) {
    fd <- central_difference(p, j, "birth_rate", NULL, fitness_of, rel = 1e-5)
    for (metric in metrics) {
      expect_equal(g$gradient[metric, paste0(j, ".birth_rate")], fd[[metric]],
                   tolerance = 1e-5, label = paste(metric, j))
    }
  }
  # offspring_i = b_i R0_i, so its own birth-rate column is R0_i plus b_i times
  # the net reproduction ratio's, and the cross column is b_i times the cross.
  for (i in 1:2) for (j in 1:2) {
    own <- g$gradient[paste0(i, ".net_reproduction_ratio"), paste0(j, ".birth_rate")]
    expected <- b[i] * own + (if (i == j) g$value[[paste0(i, ".net_reproduction_ratio")]] else 0)
    expect_equal(g$gradient[paste0(i, ".offspring_production"), paste0(j, ".birth_rate")],
                 expected, tolerance = 1e-10, label = paste(i, j))
  }
})

test_that("a variable birth rate has no column rather than a zero", {
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- 3
  p <- add_strategies(p, trait_matrix(c(0.0825, 3), c("lma", "hmat")),
                      birth_rate = list(list(x = c(0, 3), y = c(2, 2))))
  expect_true(p$strategies[[1]]$is_variable_birth_rate)
  solver <- ff16_stand_run(p)
  g <- stand_gradient(solver, metrics = "1.net_reproduction_ratio", hyperpar = NULL)
  expect_false("1.birth_rate" %in% colnames(g$gradient))
  expect_error(stand_gradient(solver, traits = "1.birth_rate", hyperpar = NULL), "Unknown trait")
})

test_that("the birth-rate column passes through a hyperparameter function as itself", {
  solver <- ff16_stand_run(ff16_stand_parameters(lifetime = 6),
                       Control(node_density_coordinate = "birth_date"))
  through <- stand_gradient(solver, metrics = "1.net_reproduction_ratio",
                            traits = c("1.lma", "1.birth_rate"),
                            hyperpar = default_hyperpar("FF16"))
  partial <- stand_gradient(solver, metrics = "1.net_reproduction_ratio",
                            traits = "1.birth_rate", hyperpar = NULL)
  expect_identical(through$gradient[, "1.birth_rate"],
                   partial$gradient[, "1.birth_rate"])
})
