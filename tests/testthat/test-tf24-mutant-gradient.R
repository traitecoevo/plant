# TF24 invaders against a recorded resident: the same replay and sweep as FF16's
# (test-ff16-mutant-gradient.R), on a model whose invaders solve their own leaf
# in the recorded field and whose sweep loads those solves.

tf24_resident <- function(lifetime = 2) {
  p <- scm_base_parameters("TF24")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = 1)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, lifetime, 0.25))
  p
}
tf24_env <- function() {
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 1)
  env
}
tf24_control <- function() {
  Control(node_density_in_birth_date = TRUE, ode_step_size_max = 0.05)
}

test_that("a TF24 invader's census gradient is the difference across replays", {
  p <- tf24_resident()
  scm <- run_scm(p, tf24_env(), tf24_control(), record_trajectory = TRUE)
  pm <- p
  pm$strategies[[1]]$pars$lma <- 0.1
  scm$run_mutant(pm)
  metrics <- c("leaf_area", "1.net_reproduction_ratio")
  g <- stand_gradient(scm, metrics = metrics, traits = "1.lma", hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  x <- pm$strategies[[1]]$pars$lma
  h <- 1e-5 * x
  at <- function(v) {
    q <- pm
    q$strategies[[1]]$pars$lma <- v
    scm$run_mutant(q)
    c(stand_census(scm), stand_fitness(scm))[metrics]
  }
  fd <- (at(x + h) - at(x - h)) / (2 * h)
  for (metric in metrics) {
    expect_lt(abs(g$gradient[metric, "1.lma"] - fd[[metric]]),
              1e-3 * abs(fd[[metric]]), label = metric)
  }
})

test_that("a TF24 invader solves its own leaf, alone or beside another", {
  p <- tf24_resident()
  scm <- run_scm(p, tf24_env(), tf24_control())
  resident <- scm$net_reproduction_ratios
  invaders <- function(lma) {
    q <- add_strategies(scm_base_parameters("TF24"), trait_matrix(lma, "lma"),
                        hyperpar = TF24_hyperpar,
                        birth_rate = rep(1, length(lma)))
    q$patch_type <- p$patch_type
    q$max_patch_lifetime <- p$max_patch_lifetime
    q$node_schedule_times <- rep(p$node_schedule_times, length(lma))
    q
  }
  scm$run_mutant(invaders(c(0.0825, 0.1)))
  together <- scm$net_reproduction_ratios
  scm$run_mutant(invaders(0.1))
  alone <- scm$net_reproduction_ratios
  # The resident as its own invader solves the leaf the resident solved.
  expect_equal(together[[1]], resident, tolerance = 1e-12)
  # And an invader's fitness does not depend on who invades beside it.
  expect_identical(together[[2]], alone)
  expect_gt(abs(log(alone) - log(resident)), 0.01)
})

test_that("a TF24 invasion of a resumed resident is differentiated", {
  p <- tf24_resident()
  burn <- run_scm(p, tf24_env(), tf24_control(), run = FALSE)
  burn$run_to(1)
  state <- export_patch_state(burn)
  scm <- run_scm(set_initial_state(p, state), tf24_env(), tf24_control(),
                 record_trajectory = TRUE)
  pm <- p
  pm$initial_time <- state$time
  pm$node_schedule_times <- state$node_schedule_times
  scm$run_mutant(pm)
  g <- stand_gradient(scm, metrics = "1.net_reproduction_ratio", traits = "1.lma",
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  scm$record_trajectory <- FALSE
  x <- pm$strategies[[1]]$pars$lma
  h <- 1e-4 * x
  at <- function(v) {
    q <- pm
    q$strategies[[1]]$pars$lma <- v
    scm$run_mutant(q)
    scm$net_reproduction_ratios
  }
  expect_equal(g$gradient[["1.net_reproduction_ratio", "1.lma"]],
               (at(x + h) - at(x - h)) / (2 * h), tolerance = 1e-3)
})
