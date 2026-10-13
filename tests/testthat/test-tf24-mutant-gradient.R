# TF24 invaders against a recorded resident: the same replay and sweep as FF16's
# (test-ff16-mutant-gradient.R), on a model whose replay also loads recorded
# leaf solves. TF24 declares no fitness, so its invaders' census is what is
# differentiated.

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
  expect_error(stand_fitness(scm), "no differentiable fitness")
  pm <- p
  pm$strategies[[1]]$pars$lma <- 0.1
  scm$run_mutant(pm)
  g <- stand_gradient(scm, metrics = "leaf_area", traits = "1.lma")
  expect_false(any(stand_gradient_refused(g)))
  x <- pm$strategies[[1]]$pars$lma
  h <- 1e-5 * x
  at <- function(v) {
    q <- pm
    q$strategies[[1]]$pars$lma <- v
    scm$run_mutant(q)
    stand_census(scm)[["leaf_area"]]
  }
  fd <- (at(x + h) - at(x - h)) / (2 * h)
  expect_lt(abs(g$gradient["leaf_area", "1.lma"] - fd), 1e-3 * abs(fd))
})

test_that("a TF24 invasion of a resumed resident is refused, not mispaired", {
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
  expect_error(stand_gradient(scm, metrics = "leaf_area", traits = "1.lma"),
               "resumed resident")
})
