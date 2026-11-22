# stand_gradient(mode = "forward") on FF16: forward tangents of the recorded
# run, one trait column per lane and up to eight lanes a walk, refereed by the
# reverse sweep of the same run. The two differentiate one function at the same
# steps, so they agree to rounding; the shape, names and refusals of the answer
# are the sweep's exactly.
#
# The stands are helper-ff16.R's.

# The forward answer `f` against the sweep's `r`: every element but the numbers
# identical, and the numbers each within `tolerance` of its metric's largest
# entry, so a column the sweep finds exactly zero is not read against itself.
expect_same_gradient <- function(f, r, tolerance = 1e-12) {
  expect_identical(names(f), names(r))
  expect_identical(dimnames(f$gradient), dimnames(r$gradient))
  expect_identical(f$value, r$value)
  expect_identical(f$refusal, r$refusal)
  expect_identical(f$control, r$control)
  scale <- pmax(apply(abs(r$gradient), 1, max), .Machine$double.xmin)
  expect_lte(max(abs(f$gradient - r$gradient) / scale), tolerance)
}

both_modes <- function(solver, ...) {
  list(forward = stand_gradient(solver, ..., mode = "forward"),
       reverse = stand_gradient(solver, ..., mode = "reverse"))
}

# Invaders at `lma` of the resident `solver`, on its node schedule.
forward_invaders <- function(solver, lma, hmat = 3) {
  p <- solver$parameters
  p$strategies <- ff16_stand_parameters(lma, hmat)$strategies
  p$node_schedule_times <- rep(p$node_schedule_times[1], length(lma))
  p
}

test_that("the forward route is the sweep on a resident stand", {
  solver <- ff16_stand_run()
  columns <- census_trait_names(solver, "FF16")
  # One column; three, in one walk at four lanes; ten, a walk of eight and one
  # of two; and every column, five walks.
  for (traits in list("1.lma", c("1.lma", "1.hmat", "1.omega"), columns[1:10],
                      NULL)) {
    g <- both_modes(solver, traits = traits, hyperpar = NULL)
    expect_false(any(stand_gradient_refused(g$forward)))
    expect_same_gradient(g$forward, g$reverse)
  }
  expect_identical(colnames(both_modes(solver, hyperpar = NULL)$forward$gradient),
                   columns)
})

test_that("the forward route is the sweep for invaders, across walks", {
  solver <- ff16_stand_run()
  solver$run_mutant(forward_invaders(solver, 0.07))
  metrics <- c("1.net_reproduction_ratio", "1.offspring_production")
  columns <- census_trait_names(solver, "FF16")
  # Nine columns: a walk of eight, then a walk of one.
  for (traits in list("1.lma", c("1.lma", "1.omega", "1.a_f1"), columns[1:9])) {
    g <- both_modes(solver, metrics = metrics, traits = traits, hyperpar = NULL)
    expect_same_gradient(g$forward, g$reverse)
  }
  # Two invaders: each one's columns are its own, species-major.
  solver$run_mutant(forward_invaders(solver, c(0.07, 0.1)))
  g <- both_modes(solver, metrics = c("1.net_reproduction_ratio",
                                      "2.net_reproduction_ratio"),
                  traits = c("2.hmat", "1.lma", "2.lma"), hyperpar = NULL)
  expect_same_gradient(g$forward, g$reverse)
  expect_identical(g$forward$gradient["1.net_reproduction_ratio", "2.lma"], 0)
})

test_that("the hyperparameter chain is applied alike in both modes", {
  solver <- ff16_stand_run()
  traits <- c("1.lma", "1.hmat")
  g <- both_modes(solver, traits = traits)
  expect_same_gradient(g$forward, g$reverse)
  # The chain reached lma: its total is not its partial.
  partial <- stand_gradient(solver, traits = traits, hyperpar = NULL)
  expect_false(isTRUE(all.equal(g$forward$gradient[, "1.lma"],
                                partial$gradient[, "1.lma"], tolerance = 1e-3)))
  # FF16's trade-off derives nine parameters, and only k_l and r_l move with
  # lma or hmat, so "auto" walks four columns forward rather than eleven.
  expect_identical(stand_gradient(solver, traits = traits), g$forward)
  # Every trait the strategy has, which the sweep takes under "auto".
  g <- both_modes(solver, metrics = "1.offspring_production")
  expect_same_gradient(g$forward, g$reverse)
})

test_that("auto takes the forward route for FF16 below the crossover only", {
  solver <- ff16_stand_run(ff16_stand_parameters(lifetime = 3))
  expect_true(use_forward_gradient(solver, "auto", forward_gradient_crossover - 1L))
  expect_false(use_forward_gradient(solver, "auto", forward_gradient_crossover))
  expect_false(use_forward_gradient(solver, "reverse", 1L))
  expect_true(use_forward_gradient(solver, "forward", 100L))
  traits <- c("1.lma", "1.hmat", "1.omega")
  expect_identical(stand_gradient(solver, traits = traits, hyperpar = NULL),
                   stand_gradient(solver, traits = traits, hyperpar = NULL,
                                  mode = "forward"))
  expect_identical(stand_gradient(solver, hyperpar = NULL),
                   stand_gradient(solver, hyperpar = NULL, mode = "reverse"))
  # A repeated column is computed once and reported where it was asked for,
  # as the sweep reports it, and so is a repeated metric.
  g <- both_modes(solver, metrics = c("1.offspring_production", "leaf_area",
                                      "1.offspring_production"),
                  traits = c("1.hmat", "1.lma", "1.hmat"), hyperpar = NULL)
  expect_same_gradient(g$forward, g$reverse)
})

test_that("the forward route runs again a run that kept no states, and leaves it", {
  solver <- run(ff16_stand_parameters(lifetime = 3), NULL,
                ff16_referee_control())$solver
  census <- stand_census(solver)
  f <- stand_gradient(solver, traits = c("1.lma", "1.k_I"), hyperpar = NULL,
                      mode = "forward")
  expect_true(solver$keep_states)
  expect_identical(stand_census(solver), census)
  r <- stand_gradient(solver, traits = c("1.lma", "1.k_I"), hyperpar = NULL,
                      mode = "reverse")
  expect_same_gradient(f, r)
  # And a second forward call gives the same numbers as the first.
  expect_identical(stand_gradient(solver, traits = c("1.lma", "1.k_I"),
                                  hyperpar = NULL, mode = "forward"), f)
})

test_that("the forward route refuses what the sweep refuses", {
  forward <- function(solver, ...) {
    stand_gradient(solver, ..., hyperpar = NULL, mode = "forward")
  }
  solver <- ff16_stand(lifetime = 2)
  expect_error(forward(solver, traits = "1.eta"), "No gradient exists")
  expect_error(forward(solver, metrics = "3.offspring_production"),
               "Unknown metric")
  solver <- ff16_stand(lifetime = 2, ctrl = Control())
  expect_error(forward(solver, traits = "1.lma"), "birth-date")
  solver <- ff16_stand(lifetime = 2,
                       ctrl = ff16_referee_control(shading_model = "ppa"))
  expect_error(forward(solver, traits = "1.lma"), "staircase")
})

test_that("the absorbing top layer sends every call to the sweep", {
  # Once a cohort is absorbed the walk's columns through the layer (hmat,
  # a_f2) came out 1-2% off the sweep, which agrees with a pinned-step
  # difference, so the walk is refused with the layer on.
  p <- add_strategies(testing_parameters("FF16"),
                      trait_matrix(c(0.0825, 3), c("lma", "hmat")), birth_rate = 2)
  p$fitness_weighting <- "none"
  p$max_patch_lifetime <- 19
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 19, 0.5))
  solver <- run(p, NULL, ff16_referee_control(absorb_epsilon = 0.05),
                record = "states")$solver
  expect_true(any(solver$patch$species[[1]]$retired))
  traits <- c("1.lma", "1.hmat")
  expect_false(use_forward_gradient(solver, "auto", 1L))
  expect_identical(stand_gradient(solver, traits = traits, hyperpar = NULL),
                   stand_gradient(solver, traits = traits, hyperpar = NULL,
                                  mode = "reverse"))
  expect_error(stand_gradient(solver, traits = traits, hyperpar = NULL,
                              mode = "forward"),
               "does not follow cohorts into the absorbing top layer")
  expect_error(census_trait_gradient_forward(solver, "FF16", 1L),
               "does not follow cohorts into the absorbing top layer")
})

test_that("mode = \"forward\" is an error on any model but FF16", {
  p <- testing_parameters("TF24")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = list(1))
  ctrl <- Control(node_density_coordinate = "birth_date")
  # Refused before anything runs, so the runner need not have run.
  tf24 <- DeterministicPatchRunner("TF24", "TF24_Env")(p, Environment("TF24"),
                                                       empty_events(), ctrl)
  expect_error(stand_gradient(tf24, traits = "1.lma", mode = "forward"),
               "available for FF16 only, not the TF24 model: its strategy keeps")
  expect_false(use_forward_gradient(tf24, "auto", 1L))
  # The C++ entry refuses it too, and never compiles the walk for it.
  expect_error(census_trait_gradient_forward(tf24, "TF24", 1L),
               "keeps the choices of an inner solve")

  p <- testing_parameters("K93")
  p$max_patch_lifetime <- 2
  p <- add_strategies(p, trait_matrix(0.059, "b_0"), birth_rate = list(1))
  k93 <- run(p, NULL, ctrl)$solver
  expect_error(stand_gradient(k93, hyperpar = NULL, mode = "forward"),
               "available for FF16 only, not the K93 model")
})
