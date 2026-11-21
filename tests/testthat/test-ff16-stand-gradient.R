# The stand gradient on FF16: the same entry points TF24 answers, refereed by a
# central difference of the SAME discretisation.
#
# The reverse sweep transposes the steps the run recorded, so it is the
# derivative of the run at those step sizes. A difference that lets the adaptive
# stepper re-choose its steps differentiates the controller too, and scatters at
# the percent level; pinning the recorded times and sizes onto the perturbed runs
# removes that.
#
# ⚠️ THE REFEREE RUNS AT A TIGHT SEED-HEIGHT ROOT TOLERANCE, set here although
# it is the default, so a looser default cannot weaken the referee; the
# gradient does not need it. The seed's height is a root-find's iterate in the
# forward run and the exact root's derivative on the differentiated path
# (implicit function theorem). A difference at `offspring_height_tol = 1e-8`
# differentiates the iterate, whose stopping point moves with the trait
# differently from the root: measured on this stand, d(mass)/d(lma) through
# the FF16 trade-off function was off by 6.8e-4 at a relative step of 1e-4 and by 1.8e-4 at
# 1e-5, while the sweep's answer did not move with the tolerance at all. At
# 1e-12 the two agree to 2.5e-8.
referee_control <- function(...) {
  Control(node_density_coordinate = "birth_date", offspring_height_tol = 1e-12,
          ...)
}

ff16_stand <- function(lifetime = 5, lma = 0.0825, birth_rate = 2,
                       ctrl = referee_control(), ...) {
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(lma, "lma"),
                      birth_rate = as.list(rep(birth_rate, length(lma))))
  run(p, NULL, ctrl, record = "states", ...)$solver
}

census_of <- function(ctrl = referee_control()) {
  function(p) stand_census(run(p, NULL, ctrl)$solver)
}

test_that("FF16 declares the three census metrics and names its columns", {
  solver <- ff16_stand(lifetime = 2)
  expect_equal(names(stand_census(solver)),
               c("leaf_area", "mass_above_ground", "area_stem"))
  g <- stand_gradient(solver, metrics = "leaf_area", hyperpar = NULL)
  # Every parameter of FF16_Pars bar the one the model says it cannot
  # differentiate, each prefixed by its species.
  pars <- names(FF16_Model()$pars)
  cannot <- names(census_undifferentiable("FF16"))
  expect_setequal(cannot, "eta")
  # ... plus the birth rate, a column of its own for a constant rate.
  expect_setequal(colnames(g$gradient),
                  paste0("1.", c(setdiff(pars, cannot), "birth_rate")))
  expect_false(any(stand_gradient_refused(g)))
  expect_error(stand_gradient(solver, traits = "1.eta", hyperpar = NULL), "No gradient exists")
})

test_that("the FF16 census is the trapezium of density times the allometry", {
  solver <- ff16_stand(lifetime = 3)
  sp <- solver$patch$species[[1]]
  st <- do.call(cbind, c(lapply(sp$nodes, function(n) n$ode_state),
                         list(sp$new_node$ode_state)))
  rownames(st) <- sp$new_node$ode_names
  birth_date <- c(vapply(sp$nodes, function(n) n$introduction_time, 0),
                  sp$new_node$introduction_time)
  pars <- solver$parameters$strategies[[1]]$pars
  eta_c <- 1 - 2 / (1 + pars$eta) + 1 / (1 + 2 * pars$eta)
  height <- st["height", ]
  area_leaf <- (height / pars$a_l1)^(1 / pars$a_l2)
  area_sapwood <- area_leaf * pars$theta
  area_bark <- pars$a_b1 * area_leaf * pars$theta
  psi <- list(
    leaf_area = area_leaf,
    mass_above_ground = area_leaf * pars$lma +
      area_bark * height * eta_c * pars$rho +
      area_sapwood * height * eta_c * pars$rho +
      st["mass_heartwood", ],
    area_stem = area_bark + area_sapwood + st["area_heartwood", ])
  # The density each node contributes, which for the newest node and the
  # boundary node includes the establishment of the interval still open.
  density <- exp(vapply(c(sp$nodes, list(sp$new_node)),
                        function(n) n$log_density, 0))
  trapezium <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
  expected <- vapply(psi, function(v) trapezium(birth_date, density * v), 0)
  expect_equal(stand_census(solver), expected, tolerance = 1e-12)
})

test_that("the FF16 stand gradient is the pinned-step central difference", {
  solver <- ff16_stand(lifetime = 5)
  metrics <- c("leaf_area", "mass_above_ground", "area_stem")
  g <- stand_gradient(solver, metrics = metrics, hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  p <- pin_steps(solver)
  # One trait from each route a parameter takes into the stand: the allometry
  # (lma, a_l1, theta), the seed's residual (omega), the light field (k_I), the
  # leaf (a_p2), turnover (k_l), mortality (a_dG1, d_I) and establishment (a_d0).
  # Each metric on its own scale: area_stem is far smaller than the other two.
  for (name in c("lma", "a_l1", "theta", "omega", "k_I", "a_p2", "k_l",
                 "a_dG1", "d_I", "a_d0")) {
    fd <- central_difference(p, 1, name, NULL, census_of())
    for (m in metrics) {
      expect_equal(g$gradient[m, paste0("1.", name)], fd[[m]],
                   tolerance = 1e-5, label = paste(m, name))
    }
  }
})

test_that("respiration, turnover and net production differentiate in every parameter", {
  # The mass cascade's remaining parameters, with coarse roots switched on so
  # their respiration and turnover terms are live (a_cr1 is zero by default).
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- 5
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = list(2))
  p$strategies[[1]]$pars$a_cr1 <- 0.25
  solver <- run(p, NULL, referee_control(), record = "states")$solver
  g <- stand_gradient(solver, metrics = "mass_above_ground", hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  q <- pin_steps(solver)
  for (name in c("rho", "a_b1", "a_r1", "a_cr1", "r_s", "r_b", "r_r", "r_cr",
                 "k_b", "k_s", "k_r", "k_cr", "a_bio", "a_y")) {
    fd <- central_difference(q, 1, name, "mass_above_ground", census_of())
    expect_equal(g$gradient["mass_above_ground", paste0("1.", name)], fd,
                 tolerance = 1e-5, label = name)
  }
})

test_that("a trait's total through the FF16 trade-off function is the difference through add_strategies", {
  solver <- ff16_stand(lifetime = 5)
  # lma sets k_l and r_l through the FF16 trade-off function, so its column carries both.
  g <- stand_gradient(solver, metrics = "mass_above_ground", traits = "1.lma")
  x <- 0.0825
  h <- 1e-4 * x
  at <- function(v) {
    p <- testing_parameters("FF16")
    p$max_patch_lifetime <- 5
    p <- add_strategies(p, trait_matrix(v, "lma"), birth_rate = list(2))
    p$ode_times <- solver$ode_times
    p$ode_step_sizes <- solver$ode_step_sizes
    stand_census(run(p, NULL, referee_control())$solver)[["mass_above_ground"]]
  }
  fd <- (at(x + h) - at(x - h)) / (2 * h)
  expect_equal(g$gradient["mass_above_ground", "1.lma"], fd, tolerance = 1e-5)
})

test_that("both crown integrals differentiate through the light and the leaf rate", {
  # Deep crown integrates the leaf rate over the crown, so a_p1 reaches the census
  # as a parameter of the integrand; mean light integrates the light alone and
  # applies a_p1 after. lma reaches both through the light field and the height
  # the crown is read at.
  for (model in c("deep-crown", "mean-light")) {
    ctrl <- referee_control(shading_model = model)
    solver <- ff16_stand(lifetime = 5, ctrl = ctrl)
    g <- stand_gradient(solver, metrics = "mass_above_ground",
                        traits = c("1.lma", "1.a_p1"), hyperpar = NULL)
    for (trait in c("lma", "a_p1")) {
      fd <- central_difference(pin_steps(solver), 1, trait, "mass_above_ground",
                               census_of(ctrl))
      expect_equal(g$gradient["mass_above_ground", paste0("1.", trait)], fd,
                   tolerance = 1e-5, label = paste(model, trait))
    }
  }
})

test_that("two FF16 species carry each their own columns", {
  solver <- ff16_stand(lifetime = 5, lma = c(0.0825, 0.2))
  g <- stand_gradient(solver, metrics = "mass_above_ground",
                      traits = c("1.lma", "2.lma"), hyperpar = NULL)
  p <- pin_steps(solver)
  for (i in 1:2) {
    fd <- central_difference(p, i, "lma", "mass_above_ground", census_of())
    expect_equal(g$gradient["mass_above_ground", paste0(i, ".lma")], fd,
                 tolerance = 1e-5, label = paste0(i, ".lma"))
  }
})

test_that("the gradient refuses what it cannot differentiate, by name", {
  # The height coordinate: the same refusal TF24 makes.
  solver <- ff16_stand(lifetime = 2, ctrl = Control())
  expect_error(stand_gradient(solver, hyperpar = NULL), "birth-date")
  # The PPA stepped light profile has no derivative to carry.
  solver <- ff16_stand(lifetime = 2, ctrl = referee_control(shading_model = "ppa"))
  expect_error(stand_gradient(solver, hyperpar = NULL), "staircase")
  # A model that declares no census.
  p <- testing_parameters("K93")
  p$max_patch_lifetime <- 2
  p <- add_strategies(p, trait_matrix(0.059, "b_0"), birth_rate = list(1))
  solver <- run(p, NULL, Control(node_density_coordinate = "birth_date"))$solver
  expect_error(stand_gradient(solver, hyperpar = NULL), "No census is declared for the K93 model")
})

test_that("FF16's growth cutoff is counted where suppressed cohorts reach it, and leaves no jump", {
  # At zero net production FF16 sets the height and fecundity rates to zero.
  # Both reach zero there anyway, so the cutoff is a kink, not a jump: the run
  # is continuous across it and the sweep stays the pinned run's derivative.
  # Heartwood keeps forming there; when its rates were cut too, a cohort
  # crossing the cutoff moved the run by a step the sweep did not carry, and at
  # k_I = 4 mass_above_ground's column missed a pinned difference at a relative
  # step of 1e-4 by 5% and area_stem's by 19%. Counted, so a gradient taken
  # where it binds says so.
  #
  # At k_I = 4 suppressed cohorts reach it; at the default k_I nothing does.
  counts <- function(k_I) {
    p <- testing_parameters("FF16")
    p$max_patch_lifetime <- 5
    p <- add_strategies(p, trait_matrix(c(0.0825, k_I), c("lma", "k_I")),
                        birth_rate = list(2))
    solver <- run(p, NULL, referee_control(), record = "states")$solver
    fired <- census_clamp_counts(solver, "FF16")
    metrics <- c("mass_above_ground", "area_stem")
    g <- stand_gradient(solver, metrics = metrics, traits = "1.a_p1",
                        hyperpar = NULL)
    expect_false(any(stand_gradient_refused(g)))
    list(fired = fired, swept = census_clamp_counts(solver, "FF16", TRUE),
         solver = solver, gradient = g$gradient, metrics = metrics)
  }
  shaded <- counts(4)
  expect_gt(sum(shaded$fired[, "production_cutoff"]), 0)
  expect_gt(sum(shaded$swept[, "production_cutoff"]), 0)
  # The two metrics that read heartwood agree with the pinned difference where
  # the cutoff binds. Compared as a ratio, because area_stem's column is about
  # 1e-6 and a tolerance on the values themselves would be absolute there. The
  # step is stated, not defaulted: a heartwood rate that jumps at the cutoff
  # shows at a relative step of 1e-4 (5% and 20% off) and hides at 1e-5.
  fd <- central_difference(pin_steps(shaded$solver), 1, "a_p1", NULL, census_of(),
                           rel = 1e-4)
  for (m in shaded$metrics) {
    expect_equal(shaded$gradient[m, "1.a_p1"] / fd[[m]], 1, tolerance = 1e-5,
                 label = paste(m, "a_p1 at k_I = 4, relative to the difference"))
  }
  open <- counts(0.5)
  expect_equal(sum(open$fired[, "production_cutoff"]), 0)
  expect_equal(sum(open$swept[, "production_cutoff"]), 0)
  # FF16 reaches no other site, the light field's floor included.
  others <- setdiff(colnames(shaded$fired), "production_cutoff")
  expect_true(all(shaded$fired[, others] == 0))
  expect_true(all(shaded$swept[, others] == 0))
})
