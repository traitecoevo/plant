# Invaders of a stand at a fixed point (#519): the resident makes one field, and
# an invader stands in it, held, for its whole life (SCM$run_mutant_held,
# fixed_point_fitness). The referee is the same cohort born into the running
# stand, which sees the field the stand actually makes.

ff16_fixed_point <- function(b = 1, t_eq = 300, absorb = 0.02, by = 1) {
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- t_eq + 1
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = b)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, t_eq + 1, by))
  ctrl <- Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb)
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), ctrl)
  scm$run_to(t_eq)
  list(scm = scm, p = p, ctrl = ctrl, t_eq = t_eq, b = b)
}

# The fraction of a strategy's seeds that establish in the field the stand
# holds at the reading year.
establishment_at <- function(fp, lma) {
  p <- add_strategies(fp$p, trait_matrix(lma, "lma"), birth_rate = 1)
  s <- p$strategies[[2]]
  s$control <- fp$ctrl
  env <- fp$scm$patch$environment
  node <- Node("FF16", "FF16_Env")(s)
  node$compute_initial_conditions(env, 1, 1)
  node$individual$establishment_probability(env)
}

# The cohorts of each strategy born at the reading year, run out `life` years in
# the stand itself, per seed arriving: each cohort's output per seed that
# establishes, times the fraction that do. A window over that year exempts them
# from absorption.
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
  c(resident = sp[[1]]$net_reproduction_ratio_by_node[sp[[1]]$node_times == fp$t_eq] *
      S_D * establishment_at(fp, fp$p$strategies[[1]]$pars$lma),
    invader = sp[[2]]$net_reproduction_ratio_by_node[[1]] * S_D *
      establishment_at(fp, lma))
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
  # TF24's soil is part of the field. Held, the patch stands in it whatever state
  # it is set to, and the soil's own rates are zero, so it does not move either.
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
  pa <- scm$patch
  n_env <- pa$environment$ode_size
  y <- pa$ode_state
  soil <- tail(y, n_env)
  expect_true(any(tail(pa$ode_rates, n_env) != 0))
  pa$hold_field()
  expect_true(pa$holds_field)
  expect_true(all(tail(pa$ode_rates, n_env) == 0))
  # A state with other soil water: the patch still stands in the held soil.
  y_dry <- y
  y_dry[length(y) - n_env + 1] <- soil[[1]] / 2
  pa$set_ode_state(y_dry, pa$time)
  expect_identical(tail(pa$ode_state, n_env), soil)
  pa$release_field()
  expect_false(pa$holds_field)
  expect_true(any(tail(pa$ode_rates, n_env) != 0))
  # Through the SCM, the invader runs in it and the SCM releases it afterwards.
  R <- fixed_point_fitness(scm, p, trait_matrix(0.0825, "lma"), lifetime = 200,
                           survival_tol = 1)
  expect_true(is.finite(R))
  expect_false(scm$patch$holds_field)
})

test_that("a held field is released once the invaders have run", {
  # A later run of the SCM builds its own field: a rerun of the invaders'
  # parameters is a fresh run of them, not another held one.
  fp <- ff16_fixed_point(t_eq = 150)
  held <- fixed_point_fitness(fp$scm, fp$p, trait_matrix(0.1, "lma"))
  fp$scm$run()
  rerun <- fp$scm$patch$species[[1]]$net_reproduction_ratio_by_node[[1]] *
    fp$p$strategies[[1]]$pars$S_D
  pm <- fp$scm$parameters
  fresh <- run_scm(pm, ctrl = fp$ctrl)
  expect_identical(rerun, fresh$patch$species[[1]]$net_reproduction_ratio_by_node[[1]] *
                     fp$p$strategies[[1]]$pars$S_D)
  expect_false(isTRUE(all.equal(rerun, held)))
})

test_that("a field held inside an interval scores as one held at its start", {
  # Between two introductions the field swings; on a grid fine enough for a
  # fixed point the swing is small, which is the check a coarse grid fails.
  fp <- ff16_fixed_point(t_eq = 150)
  at_start <- fixed_point_fitness(fp$scm, fp$p, trait_matrix(0.1, "lma"))
  p <- fp$p
  p$max_patch_lifetime <- 150.5
  mid <- run_scm(p, ctrl = fp$ctrl)
  at_mid <- fixed_point_fitness(mid, p, trait_matrix(0.1, "lma"))
  expect_equal(at_mid, at_start, tolerance = 1e-2)
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
  # Survival since establishment, not since seed: in a dense stand a seedling's
  # establishment is far below any tolerance, and a cohort cut off before it
  # reproduces must still count as alive.
  dense <- ff16_fixed_point(b = 68.867, t_eq = 150)
  expect_error(fixed_point_fitness(dense$scm, dense$p, trait_matrix(0.06, "lma"),
                                   lifetime = 40),
               "still alive")
})

## The parameters fixed_point_fitness() builds for one invader.
held_invader <- function(fp, lma, lifetime = 200) {
  pm <- add_mutant(fp$p, trait_matrix(lma, "lma"), birth_rate = list(1))
  pm$patch_type <- "fixed"
  pm$max_patch_lifetime <- lifetime
  pm$node_schedule_times <- list(0)
  pm$initial_time <- 0
  pm$ode_times <- numeric(0)
  pm$ode_step_sizes <- numeric(0)
  pm
}

test_that("a gradient refuses a run made in a held field, after it is released", {
  ## The field is released when run_mutant_held() returns, but the outputs, and
  ## any recording, were made in it.
  fp <- ff16_fixed_point(t_eq = 60, absorb = 0)
  fp$scm$run_mutant_held(held_invader(fp, 0.1))
  expect_false(fp$scm$patch$holds_field)
  expect_error(stand_gradient(fp$scm, metrics = "leaf_area", hyperpar = NULL), "held field")
  fp <- ff16_fixed_point(t_eq = 60, absorb = 0)
  fp$scm$record_trajectory <- TRUE
  fp$scm$run_mutant_held(held_invader(fp, 0.1))
  expect_error(stand_gradient(fp$scm, metrics = "leaf_area", hyperpar = NULL), "held field")
  ## A later run of the SCM is its own, and the gradient takes it.
  fp$scm$run()
  expect_true(all(is.finite(stand_gradient(fp$scm, metrics = "leaf_area",
                                           hyperpar = NULL)$gradient)))
})

test_that("a held field is released when the invaders' run fails", {
  fp <- ff16_fixed_point(t_eq = 60)
  ## With the layer on, a zero birth rate is refused at the introduction, which
  ## is inside the invaders' run.
  pm <- held_invader(fp, 0.1)
  pm$strategies[[1]]$birth_rate_y <- 0
  expect_error(fp$scm$run_mutant_held(pm), "birth rate zero")
  expect_false(fp$scm$patch$holds_field)
  ## And an invader is checked as a resident is: with no asymptote to absorb at,
  ## it is refused before it runs.
  pm <- held_invader(fp, 0.1)
  pm$strategies[[1]]$pars$a_f2 <- -1
  expect_error(fp$scm$run_mutant_held(pm), "a_f2")
  expect_false(fp$scm$patch$holds_field)
})

test_that("a coarse grid makes the fixed point a periodic orbit", {
  ## Between two introductions the field swings. At the resident's equilibrium
  ## birth rate and a cohort a year the swing is percent-level, so a field held
  ## at an introduction is not the stand's; a cohort a quarter year shrinks it
  ## a hundredfold. Read with the resident's own strategy.
  b_star <- 68.867
  swing <- function(by) {
    fp <- ff16_fixed_point(b = b_star, t_eq = 150, by = by)
    at_start <- fixed_point_fitness(fp$scm, fp$p, trait_matrix(0.0825, "lma"))
    p <- fp$p
    p$max_patch_lifetime <- 150 + by / 2
    tt <- p$node_schedule_times[[1]]
    p$node_schedule_times <- list(tt[tt < p$max_patch_lifetime])
    mid <- run_scm(p, ctrl = fp$ctrl)
    at_mid <- fixed_point_fitness(mid, p, trait_matrix(0.0825, "lma"))
    abs(at_mid / at_start - 1)
  }
  coarse <- swing(1)
  fine <- swing(0.25)
  expect_gt(coarse, 5e-3)
  expect_lt(fine, 1e-3)
  expect_lt(fine, coarse / 10)
})
