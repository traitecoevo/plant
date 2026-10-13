# Invaders of a stand at a fixed point (#519): the resident makes one field, and
# an invader stands in it, held, for its whole life (SCM$run_mutant_held,
# fixed_point_fitness). The referee is the same cohort born into the running
# stand, which sees the field the stand actually makes.

ff16_fixed_point <- function(b = 1, t_eq = 300, absorb = 0.02) {
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- t_eq + 1
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = b)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, t_eq + 1, 1))
  ctrl <- Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb)
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), ctrl)
  scm$run_to(t_eq)
  list(scm = scm, p = p, ctrl = ctrl, t_eq = t_eq, b = b)
}

# The cohorts of each strategy born at the reading year, run out `life` years in
# the stand itself. A window over that year exempts them from absorption.
in_stand_fitness <- function(fp, lma, life = 400) {
  p <- fp$p
  p$max_patch_lifetime <- fp$t_eq + life
  p$node_schedule_times <- list(node_schedule_times_periodic(0, fp$t_eq + life, 1))
  p <- add_strategies(p, trait_matrix(lma, "lma"),
                      birth_rate = fp$b * 1e-6)
  p$node_schedule_times[[2]] <- fp$t_eq
  p <- fitness_window(p, fp$t_eq, fp$t_eq + 1)
  scm <- run_scm(p, ctrl = fp$ctrl)
  sp <- scm$patch$species
  S_D <- p$strategies[[1]]$pars$S_D
  c(resident = sp[[1]]$net_reproduction_ratio_by_node[sp[[1]]$node_times == fp$t_eq] * S_D,
    invader = sp[[2]]$net_reproduction_ratio_by_node[[1]] * S_D)
}

test_that("an invader held in the fixed point's field scores its cohort in the stand", {
  fp <- ff16_fixed_point()
  lma <- c(0.0825, 0.06, 0.1)
  held <- fixed_point_fitness(fp$scm, fp$p, trait_matrix(lma, "lma"))
  expect_true(all(is.finite(held)) && all(held > 0))
  for (k in seq_along(lma)) {
    ref <- in_stand_fitness(ff16_fixed_point(), lma[k])
    expect_equal(held[[k]], ref[["invader"]], tolerance = 2e-3)
  }
  # The resident's own strategy scores the resident's own cohort.
  expect_equal(held[[1]], in_stand_fitness(ff16_fixed_point(), 0.0825)[["resident"]],
               tolerance = 2e-3)
  # A thicker-leaved invader does better in this closed canopy, a thinner one
  # worse: the field is the resident's, not an empty patch's.
  expect_gt(held[[3]], held[[1]])
  expect_lt(held[[2]], held[[1]])
})

test_that("every invader stands in the field taken at the first call", {
  fp <- ff16_fixed_point(t_eq = 150)
  a <- fixed_point_fitness(fp$scm, fp$p, trait_matrix(0.1, "lma"))
  # Run again alone, the same strategy takes the same steps in the same field.
  expect_identical(fixed_point_fitness(fp$scm, fp$p, trait_matrix(0.1, "lma")), a)
  # Beside another it shares their step sequence, so it agrees to the solver's
  # tolerance rather than bit for bit.
  b <- fixed_point_fitness(fp$scm, fp$p, trait_matrix(c(0.06, 0.1), "lma"))
  expect_equal(b[[2]], a, tolerance = 1e-3)
})

test_that("a held field holds the environment's own state", {
  # TF24's soil is part of the field; held, it does not move however long the
  # invader runs.
  p <- scm_base_parameters("TF24")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- 6
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = 1)
  p$node_schedule_times <- list(0:5)
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 1)
  ctrl <- Control(node_density_in_birth_date = TRUE, ode_step_size_max = 0.05)
  scm <- SCM("TF24", "TF24_Env")(p, env, empty_events(), ctrl)
  scm$run_to(5)
  n_env <- scm$patch$environment$ode_size
  soil <- tail(scm$patch$ode_state, n_env)
  R <- fixed_point_fitness(scm, p, trait_matrix(0.0825, "lma"), lifetime = 200,
                           survival_tol = 1)
  expect_true(is.finite(R))
  expect_identical(tail(scm$patch$ode_state, n_env), soil)
})

test_that("holding a field needs a resident to take it from", {
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(),
                                 Control(node_density_in_birth_date = TRUE))
  expect_error(scm$run_mutant_held(p), "fixed point")
})

test_that("an invader that outlives its run is refused", {
  fp <- ff16_fixed_point(t_eq = 100)
  expect_error(fixed_point_fitness(fp$scm, fp$p, trait_matrix(0.0825, "lma"),
                                   lifetime = 20),
               "still alive")
})
