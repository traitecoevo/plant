# The stand gradient on FF16: the same entry points TF24 answers, refereed by a
# central difference of the SAME discretisation.
#
# The reverse sweep transposes the steps the run recorded, so it is the
# derivative of the run at those step sizes. A difference that lets the adaptive
# stepper re-choose its steps differentiates the controller too, and scatters at
# the percent level; pinning the recorded times and sizes onto the perturbed runs
# removes that.
#
# ⚠️ THE REFEREE RUNS AT A TIGHT SEED-HEIGHT ROOT TOLERANCE, and the gradient
# does not need it. The seed's height is a bisection iterate in the forward run
# and the exact root's derivative on the differentiated path (implicit function
# theorem). A difference at the default `offspring_production_tol = 1e-8`
# differentiates the iterate, whose stopping point moves with the trait
# differently from the root: measured on this stand, d(mass)/d(lma) through
# FF16_hyperpar was off by 6.8e-4 at a relative step of 1e-4 and by 1.8e-4 at
# 1e-5, while the sweep's answer did not move with the tolerance at all. At
# 1e-12 the two agree to 2.5e-8.
referee_control <- function(...) {
  Control(node_density_in_birth_date = TRUE, offspring_production_tol = 1e-12,
          ...)
}

ff16_stand <- function(lifetime = 5, lma = 0.0825, birth_rate = 2,
                       ctrl = referee_control(), ...) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(lma, "lma"),
                      birth_rate = as.list(rep(birth_rate, length(lma))))
  run_scm(p, NULL, ctrl, record_trajectory = TRUE, ...)
}

# The perturbed run stepped exactly where the recorded one stepped.
pinned <- function(scm) {
  p <- scm$parameters
  p$ode_times <- scm$ode_times
  p$ode_step_sizes <- scm$ode_step_sizes
  p
}

central_difference <- function(p, species, name, metric, rel = 1e-4,
                               ctrl = referee_control()) {
  x <- p$strategies[[species]]$pars[[name]]
  h <- rel * abs(x)
  at <- function(v) {
    q <- p
    q$strategies[[species]]$pars[[name]] <- v
    stand_census(run_scm(q, NULL, ctrl))[[metric]]
  }
  (at(x + h) - at(x - h)) / (2 * h)
}

test_that("FF16 declares the three census metrics and names its columns", {
  scm <- ff16_stand(lifetime = 2)
  expect_equal(names(stand_census(scm)),
               c("leaf_area", "mass_above_ground", "area_stem"))
  g <- stand_gradient(scm, metrics = "leaf_area", hyperpar = NULL)
  # Every parameter of FF16_Pars bar the one the model says it cannot
  # differentiate, each prefixed by its species.
  pars <- names(FF16_Strategy()$pars)
  cannot <- names(census_undifferentiable_ff16())
  expect_setequal(cannot, "eta")
  # ... plus the birth rate, a column of its own for a constant rate.
  expect_setequal(colnames(g$gradient),
                  paste0("1.", c(setdiff(pars, cannot), "birth_rate")))
  expect_false(any(stand_gradient_refused(g)))
  expect_error(stand_gradient(scm, traits = "1.eta", hyperpar = NULL), "No gradient exists")
})

test_that("the FF16 census is density times the allometry by establishment weight", {
  scm <- ff16_stand(lifetime = 3)
  sp <- scm$patch$species[[1]]
  st <- do.call(cbind, c(lapply(sp$nodes, function(n) n$ode_state),
                         list(sp$new_node$ode_state)))
  rownames(st) <- sp$new_node$ode_names
  pars <- scm$parameters$strategies[[1]]$pars
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
  density <- exp(vapply(c(sp$nodes, list(sp$new_node)),
                        function(n) n$log_density, 0))
  weight <- ladder_establishment_weights(sp)
  expected <- vapply(psi, function(v) sum(weight * density * v), 0)
  expect_equal(stand_census(scm), expected, tolerance = 1e-12)
})

test_that("the FF16 stand gradient is the pinned-step central difference", {
  scm <- ff16_stand(lifetime = 5)
  g <- stand_gradient(scm, metrics = "mass_above_ground", hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  p <- pinned(scm)
  # One trait from each route a parameter takes into the stand: the allometry
  # (lma, a_l1, theta), the seed's residual (omega), the light field (k_I), the
  # leaf (a_p2), turnover (k_l), mortality (a_dG1, d_I) and establishment (a_d0).
  for (name in c("lma", "a_l1", "theta", "omega", "k_I", "a_p2", "k_l",
                 "a_dG1", "d_I", "a_d0")) {
    fd <- central_difference(p, 1, name, "mass_above_ground")
    expect_equal(g$gradient["mass_above_ground", paste0("1.", name)], fd,
                 tolerance = 1e-5, label = name)
  }
})

test_that("a trait's total through FF16_hyperpar is the difference through add_strategies", {
  scm <- ff16_stand(lifetime = 5)
  # lma sets k_l and r_l through FF16_hyperpar, so its column carries both.
  g <- stand_gradient(scm, metrics = "mass_above_ground", traits = "1.lma")
  x <- 0.0825
  h <- 1e-4 * x
  at <- function(v) {
    p <- scm_base_parameters("FF16")
    p$max_patch_lifetime <- 5
    p <- add_strategies(p, trait_matrix(v, "lma"), birth_rate = list(2))
    p$ode_times <- scm$ode_times
    p$ode_step_sizes <- scm$ode_step_sizes
    stand_census(run_scm(p, NULL, referee_control()))[["mass_above_ground"]]
  }
  fd <- (at(x + h) - at(x - h)) / (2 * h)
  expect_equal(g$gradient["mass_above_ground", "1.lma"], fd, tolerance = 1e-5)
})

test_that("two FF16 species carry each their own columns", {
  scm <- ff16_stand(lifetime = 5, lma = c(0.0825, 0.2))
  g <- stand_gradient(scm, metrics = "mass_above_ground",
                      traits = c("1.lma", "2.lma"), hyperpar = NULL)
  p <- pinned(scm)
  for (i in 1:2) {
    fd <- central_difference(p, i, "lma", "mass_above_ground")
    expect_equal(g$gradient["mass_above_ground", paste0(i, ".lma")], fd,
                 tolerance = 1e-5, label = paste0(i, ".lma"))
  }
})

test_that("the gradient refuses what it cannot differentiate, by name", {
  # The height coordinate: the same refusal TF24 makes.
  scm <- ff16_stand(lifetime = 2, ctrl = Control())
  expect_error(stand_gradient(scm, hyperpar = NULL), "birth-date")
  # The PPA stepped light profile has no derivative to carry.
  scm <- ff16_stand(lifetime = 2, ctrl = referee_control(shading_model = "ppa"))
  expect_error(stand_gradient(scm, hyperpar = NULL), "staircase")
  # A model that declares no census.
  p <- scm_base_parameters("K93")
  p$max_patch_lifetime <- 2
  p <- add_strategies(p, trait_matrix(0.059, "b_0"), birth_rate = list(1))
  scm <- run_scm(p, NULL, Control(node_density_in_birth_date = TRUE))
  expect_error(stand_gradient(scm, hyperpar = NULL), "No census is declared for the K93 model")
})
