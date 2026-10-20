# The rare-mutant derivative on FF16: an invader's fitness differentiated by
# its own traits with the canopy held at the recorded resident (plant#650, E2).
#
# run_mutant() runs the invaders on the resident's steps, reading the field the
# resident made at every rate evaluation, in the run and in its sweep alike: the
# field is exogenous to the invader and carries no derivative. The referee is
# the same invasion at perturbed traits, which steps where the resident stepped.

mutant_control <- function(...) {
  Control(node_density_in_birth_date = TRUE, offspring_production_tol = 1e-12,
          ...)
}

resident_parameters <- function(lma = 0.0825, hmat = 3, lifetime = 8,
                                birth_rate = 2) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  m <- trait_matrix(c(lma, rep(hmat, length(lma))), c("lma", "hmat"))
  add_strategies(p, m, birth_rate = as.list(rep(birth_rate, length(lma))))
}

# The invaders' parameters: the resident's, with the strategies replaced. They
# carry the resident's node schedule, which is what the replay introduces by.
mutant_parameters <- function(scm, lma, hmat = 3) {
  p <- scm$parameters
  m <- trait_matrix(c(lma, rep(hmat, length(lma))), c("lma", "hmat"))
  q <- add_strategies(scm_base_parameters("FF16"), m,
                      birth_rate = as.list(rep(2, length(lma))))
  p$strategies <- q$strategies
  p$node_schedule_times <- rep(p$node_schedule_times[1], length(lma))
  p
}

mutant_fitness <- function(scm, lma) {
  scm$run_mutant(mutant_parameters(scm, lma))
  stand_fitness(scm)
}

test_that("a mutant's fitness gradient is the pinned difference across replays", {
  scm <- run_scm(resident_parameters(), NULL, mutant_control(),
                 record_trajectory = TRUE)
  resident <- stand_gradient(scm, metrics = "1.net_reproduction_ratio",
                             traits = "1.lma")

  # One invader at the resident's own trait, replayed with its states kept.
  scm$run_mutant(mutant_parameters(scm, 0.0825))
  g <- stand_gradient(scm, metrics = "1.net_reproduction_ratio",
                      traits = c("1.lma", "1.omega", "1.a_f1", "1.k_I"))
  expect_false(any(stand_gradient_refused(g)))
  # The resident as its own invader has the resident's fitness.
  expect_equal(g$value[["1.net_reproduction_ratio"]],
               resident$value[["1.net_reproduction_ratio"]], tolerance = 1e-10)

  # Each column against a central difference of the replay itself. The
  # recording stands; only the invader's parameter moves.
  for (name in c("lma", "omega", "a_f1", "k_I")) {
    x <- scm$parameters$strategies[[1]]$pars[[name]]
    h <- 1e-5 * abs(x)
    at <- function(v) {
      p <- mutant_parameters(scm, 0.0825)
      p$strategies[[1]]$pars[[name]] <- v
      scm$run_mutant(p)
      stand_fitness(scm)[["1.net_reproduction_ratio"]]
    }
    fd <- (at(x + h) - at(x - h)) / (2 * h)
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
  scm <- run_scm(resident_parameters(), NULL, mutant_control(),
                 record_trajectory = TRUE)
  scm$run_mutant(mutant_parameters(scm, c(0.07, 0.1)))
  g <- stand_gradient(scm, metrics = c("1.net_reproduction_ratio",
                                       "2.net_reproduction_ratio"),
                      traits = c("1.lma", "2.lma"))
  expect_false(any(stand_gradient_refused(g)))
  for (i in 1:2) {
    x <- scm$parameters$strategies[[i]]$pars$lma
    h <- 1e-5 * x
    at <- function(v) {
      p <- mutant_parameters(scm, c(0.07, 0.1))
      p$strategies[[i]]$pars$lma <- v
      scm$run_mutant(p)
      stand_fitness(scm)[[paste0(i, ".net_reproduction_ratio")]]
    }
    fd <- (at(x + h) - at(x - h)) / (2 * h)
    expect_equal(g$gradient[paste0(i, ".net_reproduction_ratio"),
                            paste0(i, ".lma")], fd,
                 tolerance = 1e-5, label = paste0(i, ".lma"))
  }
})

test_that("a sweep of an invasion that kept no states runs it again, as an invasion", {
  scm <- run_scm(resident_parameters(lifetime = 3), NULL, mutant_control())
  scm$run_mutant(mutant_parameters(scm, 0.0825))
  again <- stand_gradient(scm, metrics = "1.net_reproduction_ratio")
  scm$record_trajectory <- TRUE
  scm$run_mutant(mutant_parameters(scm, 0.0825))
  kept <- stand_gradient(scm, metrics = "1.net_reproduction_ratio")
  expect_identical(again$gradient, kept$gradient)
})

# d(metric)/d(name) of invader `i` by central difference across replays, each
# replay of `p` with that one parameter moved.
replay_difference <- function(scm, p, i, name, metric, rel = 1e-5) {
  x <- p$strategies[[i]]$pars[[name]]
  h <- rel * abs(x)
  at <- function(v) {
    q <- p
    q$strategies[[i]]$pars[[name]] <- v
    scm$run_mutant(q)
    stand_fitness(scm)[[metric]]
  }
  (at(x + h) - at(x - h)) / (2 * h)
}

test_that("invaders do not shade each other, and an invader's birth rate reaches only its output", {
  scm <- run_scm(resident_parameters(), NULL, mutant_control(),
                 record_trajectory = TRUE)
  scm$run_mutant(mutant_parameters(scm, c(0.07, 0.1)))
  g <- stand_gradient(scm, metrics = c("1.net_reproduction_ratio",
                                       "1.offspring_production"),
                      traits = c("2.lma", "2.k_I", "1.birth_rate"))
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

test_that("an invader introduced after the replay starts has a gradient", {
  # The invader's patch is empty, and FF16's environment carries no state, so
  # the replay steps with no state until its first cohort; the sweep passes
  # those steps over (odelia 0.7.3).
  scm <- run_scm(resident_parameters(), NULL, mutant_control(),
                 record_trajectory = TRUE)
  p <- mutant_parameters(scm, 0.09)
  tt <- p$node_schedule_times[[1]]
  p$node_schedule_times <- list(tt[tt >= 1])
  scm$run_mutant(p)
  g <- stand_gradient(scm, metrics = "1.net_reproduction_ratio",
                      traits = c("1.lma", "1.a_d0"))
  expect_false(any(stand_gradient_refused(g)))
  for (name in c("lma", "a_d0")) {
    fd <- replay_difference(scm, p, 1, name, "1.net_reproduction_ratio")
    expect_lt(abs(g$gradient["1.net_reproduction_ratio", paste0("1.", name)] - fd),
              1e-5 * abs(fd), label = name)
  }
})

test_that("an invader of a resumed resident is born in the resident's canopy at the resume instant", {
  # The resident resumes at 6 and its first cohort after the resume is born at
  # 6 itself, in the field of the state the resident resumed from. Establishment
  # reads it, so a_d0 is the column that shows the wrong canopy (it was 26 % off).
  p <- resident_parameters(lifetime = 12)
  p$patch_type <- "fixed"
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 12, 0.5))
  burn <- run_scm(p, NULL, mutant_control(), run = FALSE)
  burn$run_to(6)
  state <- export_patch_state(burn)
  p2 <- set_initial_state(p, state)
  scm <- run_scm(p2, NULL, mutant_control(), record_trajectory = TRUE)
  pm <- mutant_parameters(scm, 0.09)
  pm$initial_state <- numeric(0)
  pm$n_initial_cohorts <- integer(0)
  pm$initial_node_times <- numeric(0)
  pm$initial_patch_density <- numeric(0)
  pm$initial_pr_patch_survival <- numeric(0)
  pm$initial_time <- state$time
  pm$node_schedule_times <- state$node_schedule_times
  expect_equal(pm$node_schedule_times[[1]][[1]], state$time)
  scm$run_mutant(pm)
  g <- stand_gradient(scm, metrics = "1.net_reproduction_ratio",
                      traits = c("1.lma", "1.a_d0", "1.omega"))
  expect_false(any(stand_gradient_refused(g)))
  for (name in c("lma", "a_d0", "omega")) {
    fd <- replay_difference(scm, pm, 1, name, "1.net_reproduction_ratio")
    expect_lt(abs(g$gradient["1.net_reproduction_ratio", paste0("1.", name)] - fd),
              1e-5 * abs(fd), label = name)
  }
})

test_that("an invader's gradient through the run's harvest is the difference across invasions", {
  p <- resident_parameters()
  scm <- run_scm(p, NULL, mutant_control(), record_trajectory = TRUE,
                 events = events(events_default(p),
                                 harvest(time = 4, fraction = 0.5)))
  pm <- mutant_parameters(scm, 0.09)
  scm$run_mutant(pm)
  g <- stand_gradient(scm, metrics = "1.net_reproduction_ratio",
                      traits = c("1.lma", "1.a_d0"))
  expect_false(any(stand_gradient_refused(g)))
  for (name in c("lma", "a_d0")) {
    fd <- replay_difference(scm, pm, 1, name, "1.net_reproduction_ratio")
    expect_lt(abs(g$gradient["1.net_reproduction_ratio", paste0("1.", name)] - fd),
              1e-5 * abs(fd), label = name)
  }
})
