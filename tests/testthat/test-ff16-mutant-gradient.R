# The rare-mutant derivative on FF16: an invader's fitness differentiated by
# its own traits with the canopy held at the recorded resident (plant#650, E2).
#
# run_mutant() replays the resident's recorded program with the invaders in the
# patch, loading the resident's field at every rate evaluation. Asked to keep
# its states, that replay records the field it loaded in every row, and a sweep
# of the record loads it the same way: the field is exogenous to the invader
# and carries no derivative. The referee is the same replay at perturbed
# traits, which steps where the resident stepped by construction.

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

test_that("a sweep of an invasion that kept no states is refused, not re-run", {
  scm <- run_scm(resident_parameters(lifetime = 3), NULL, mutant_control())
  scm$run_mutant(mutant_parameters(scm, 0.0825))
  expect_error(stand_gradient(scm, metrics = "1.net_reproduction_ratio"),
               "run_mutant\\(\\) kept no states")
})
