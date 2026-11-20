# The rare-mutant derivative on FF16: an invader's fitness differentiated by
# its own traits in the residents' recorded field.
#
# run_mutant() runs the invaders on the resident's steps, reading the field the
# resident made at every rate evaluation, in the run and in its sweep alike: the
# field is exogenous to the invader and carries no derivative. The referee is
# the same invasion at perturbed traits, which steps where the resident stepped.

# The invaders' parameters: the resident's, with the strategies replaced. They
# carry the resident's node schedule, which is what the pinned run introduces by.
mutant_parameters <- function(solver, lma, hmat = 3) {
  p <- solver$parameters
  p$strategies <- ff16_stand_parameters(lma, hmat)$strategies
  p$node_schedule_times <- rep(p$node_schedule_times[1], length(lma))
  p
}

# The invaders' fitness when `p` is run as the mutant set in the resident `solver`,
# which pins each run to the resident's steps.
mutant_fitness_of <- function(solver) {
  function(p) {
    solver$run_mutant(p)
    stand_fitness(solver)
  }
}

test_that("a mutant's fitness gradient is the difference across pinned runs", {
  solver <- ff16_stand_run()
  resident <- stand_gradient(solver, metrics = "1.net_reproduction_ratio",
                             traits = "1.lma", hyperpar = NULL)

  # One invader at the resident's own trait, rerun with its states kept.
  solver$run_mutant(mutant_parameters(solver, 0.0825))
  g <- stand_gradient(solver, metrics = "1.net_reproduction_ratio",
                      traits = c("1.lma", "1.omega", "1.a_f1", "1.k_I"),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  # The resident as its own invader has the resident's fitness.
  expect_equal(g$value[["1.net_reproduction_ratio"]],
               resident$value[["1.net_reproduction_ratio"]], tolerance = 1e-10)

  # Each column against a central difference of the pinned run itself. The
  # recording stands; only the invader's parameter moves.
  pm <- mutant_parameters(solver, 0.0825)
  for (name in c("lma", "omega", "a_f1", "k_I")) {
    fd <- central_difference(pm, 1, name, "1.net_reproduction_ratio",
                             mutant_fitness_of(solver), rel = 1e-5)
    expect_equal(g$gradient["1.net_reproduction_ratio", paste0("1.", name)], fd,
                 tolerance = 1e-5, label = name)
  }

  # And it is not the resident-feedback derivative: an invader does not shade
  # itself back, so the two differ.
  expect_false(isTRUE(all.equal(g$gradient["1.net_reproduction_ratio", "1.lma"],
                                resident$gradient["1.net_reproduction_ratio", "1.lma"],
                                tolerance = 1e-3)))
})

test_that("several invaders against one resident each carry their own columns", {
  solver <- ff16_stand_run()
  solver$run_mutant(mutant_parameters(solver, c(0.07, 0.1)))
  g <- stand_gradient(solver, metrics = c("1.net_reproduction_ratio",
                                          "2.net_reproduction_ratio"),
                      traits = c("1.lma", "2.lma"), hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  pm <- mutant_parameters(solver, c(0.07, 0.1))
  for (i in 1:2) {
    fd <- central_difference(pm, i, "lma", paste0(i, ".net_reproduction_ratio"),
                             mutant_fitness_of(solver), rel = 1e-5)
    expect_equal(g$gradient[paste0(i, ".net_reproduction_ratio"),
                            paste0(i, ".lma")], fd,
                 tolerance = 1e-5, label = paste0(i, ".lma"))
  }
})

test_that("a sweep of an invasion that kept no states runs it again, as an invasion", {
  solver <- run(ff16_stand_parameters(lifetime = 3), NULL, ff16_referee_control())$solver
  solver$run_mutant(mutant_parameters(solver, 0.0825))
  again <- stand_gradient(solver, metrics = "1.net_reproduction_ratio",
                          hyperpar = NULL)
  solver$keep_states <- TRUE
  solver$run_mutant(mutant_parameters(solver, 0.0825))
  kept <- stand_gradient(solver, metrics = "1.net_reproduction_ratio",
                         hyperpar = NULL)
  expect_identical(again$gradient, kept$gradient)
})

test_that("invaders do not shade each other, and an invader's birth rate reaches only its output", {
  solver <- ff16_stand_run()
  solver$run_mutant(mutant_parameters(solver, c(0.07, 0.1)))
  g <- stand_gradient(solver, metrics = c("1.net_reproduction_ratio",
                                          "1.offspring_production"),
                      traits = c("2.lma", "2.k_I", "1.birth_rate"),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  # Each stands in the resident's field, which neither moves.
  expect_identical(g$gradient["1.net_reproduction_ratio", c("2.lma", "2.k_I")],
                   c("2.lma" = 0, "2.k_I" = 0))
  # A rare invader's birth rate moves no density it sees, so its ratio does not
  # depend on it, and its offspring production is the ratio times that rate.
  expect_identical(g$gradient["1.net_reproduction_ratio", "1.birth_rate"], 0)
  expect_equal(g$gradient["1.offspring_production", "1.birth_rate"],
               g$value[["1.net_reproduction_ratio"]], tolerance = 1e-12)
})

test_that("an invader introduced after the pinned run starts has a gradient", {
  # The invader's patch is empty, and FF16's environment carries no state, so
  # the pinned run steps with no state until its first cohort; the sweep passes
  # those steps over.
  solver <- ff16_stand_run()
  p <- mutant_parameters(solver, 0.09)
  tt <- p$node_schedule_times[[1]]
  p$node_schedule_times <- list(tt[tt >= 1])
  solver$run_mutant(p)
  g <- stand_gradient(solver, metrics = "1.net_reproduction_ratio",
                      traits = c("1.lma", "1.a_d0"), hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  for (name in c("lma", "a_d0")) {
    fd <- central_difference(p, 1, name, "1.net_reproduction_ratio",
                             mutant_fitness_of(solver), rel = 1e-5)
    expect_lt(abs(g$gradient["1.net_reproduction_ratio", paste0("1.", name)] - fd),
              1e-5 * abs(fd), label = name)
  }
})

test_that("an invader of a resumed resident is born in the resident's canopy at the resume instant", {
  # The resident resumes at 6 and its first cohort after the resume is born at
  # 6 itself, in the field of the state the resident resumed from. Establishment
  # reads it, so a_d0 is the column that shows the wrong canopy.
  p <- ff16_stand_parameters(lifetime = 12)
  p$fitness_weighting <- "none"
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 12, 0.5))
  burn <- run(p, NULL, ff16_referee_control(), until = 0)$solver
  burn$run_to(6)
  state <- export_patch_state(burn)
  p2 <- set_initial_state(p, state)
  solver <- run(p2, NULL, ff16_referee_control(), record = "states")$solver
  pm <- mutant_parameters(solver, 0.09)
  pm$initial_state <- numeric(0)
  pm$n_initial_cohorts <- integer(0)
  pm$initial_node_times <- numeric(0)
  pm$initial_patch_density <- numeric(0)
  pm$initial_pr_patch_survival <- numeric(0)
  pm$initial_time <- state$time
  pm$node_schedule_times <- state$node_schedule_times
  expect_equal(pm$node_schedule_times[[1]][[1]], state$time)
  solver$run_mutant(pm)
  g <- stand_gradient(solver, metrics = "1.net_reproduction_ratio",
                      traits = c("1.lma", "1.a_d0", "1.omega"), hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  for (name in c("lma", "a_d0", "omega")) {
    fd <- central_difference(pm, 1, name, "1.net_reproduction_ratio",
                             mutant_fitness_of(solver), rel = 1e-5)
    expect_lt(abs(g$gradient["1.net_reproduction_ratio", paste0("1.", name)] - fd),
              1e-5 * abs(fd), label = name)
  }
})

test_that("an invader's gradient through the run's harvest is the difference across invasions", {
  p <- ff16_stand_parameters()
  solver <- run(p, NULL, ff16_referee_control(),
                events = events(events_default(p),
                                harvest(time = 4, fraction = 0.5)), record = "states")$solver
  pm <- mutant_parameters(solver, 0.09)
  solver$run_mutant(pm)
  g <- stand_gradient(solver, metrics = "1.net_reproduction_ratio",
                      traits = c("1.lma", "1.a_d0"), hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  for (name in c("lma", "a_d0")) {
    fd <- central_difference(pm, 1, name, "1.net_reproduction_ratio",
                             mutant_fitness_of(solver), rel = 1e-5)
    expect_lt(abs(g$gradient["1.net_reproduction_ratio", paste0("1.", name)] - fd),
              1e-5 * abs(fd), label = name)
  }
})
