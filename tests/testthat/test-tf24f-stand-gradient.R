# TF24f's stand gradient: TF24's census and fitness, with the collar's tracked
# position carried through the sweep as an ODE state. Refereed against central
# differences of runs pinned to the recorded steps, as TF24's and FF16's are
# (test-tf24-fitness-gradient.R, test-ff16-stand-gradient.R).
#
# Two stands. On the one test-model-tf24f.R runs (five years, hmat 5, a birth
# rate of 20, the default soil) the position stays near the middle of the
# feasible interval. Under dry pulses, at a gain high enough to reach it, the
# position is held at the interval's dry end through each dry spell, and the
# interval itself moves as the soil dries, so the collar's rows there are the
# end's rows.

tf24f_gradient_parameters <- function(model = "TF24f", k_acclim = 1,
                                      lifetime = 5, lma = 0.0825) {
  p <- testing_parameters(model)
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(c(lma, 5), c("lma", "hmat")),
                      hyperpar = default_hyperpar(model), birth_rate = 20)
  if (model == "TF24f") {
    s <- p$strategies[[1]]
    s$k_acclim <- k_acclim
    p$strategies[[1]] <- s
  }
  p
}

tf24f_dry_pulses <- function(lifetime) {
  env <- Environment("TF24f")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.2, 5))
  x <- seq(0, lifetime, length.out = 400)
  env$extrinsic_drivers_set_variable("rainfall", x = x,
                                     y = 0.3 * pmax(sin(2 * pi * x), 0))
  env
}

# The gradient of `metrics` and its central difference over runs pinned to the
# recorded steps, as relative errors: one row per metric, one column per trait.
tf24f_gradient_errors <- function(p, env, metrics, traits, rel = 1e-5) {
  ctrl <- Control()
  solver <- run(p, env(), ctrl, record = "states")$solver
  g <- stand_gradient(solver, metrics = metrics, traits = paste0("1.", traits),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))

  pinned <- solver$parameters
  pinned$ode_times <- solver$ode_times
  pinned$ode_step_sizes <- solver$ode_step_sizes
  pinned$ode_soil_input_slopes <- solver$ode_soil_input_slopes
  pinned$ode_soil_substep_ends <- solver$ode_soil_substep_ends
  metrics_of <- function(q) {
    s <- run(q, env(), ctrl)$solver
    c(stand_census(s), stand_fitness(s))[metrics]
  }
  expect_equal(metrics_of(pinned), g$value[metrics], tolerance = 1e-12)
  vapply(traits, function(name) {
    set <- function(v) {
      q <- pinned
      if (name == "birth_rate") {
        q$strategies[[1]]$birth_rate_y <- v
      } else {
        q$strategies[[1]]$pars[[name]] <- v
      }
      q
    }
    x <- if (name == "birth_rate") pinned$strategies[[1]]$birth_rate_y
         else pinned$strategies[[1]]$pars[[name]]
    h <- rel * abs(x)
    fd <- (metrics_of(set(x + h)) - metrics_of(set(x - h))) / (2 * h)
    abs(g$gradient[metrics, paste0("1.", name)] - fd) / abs(fd)
  }, numeric(length(metrics)))
}

test_that("TF24f declares TF24's census and fitness", {
  solver <- run(tf24f_gradient_parameters(), Environment("TF24f"), Control(),
                record = "states")$solver
  expect_equal(stand_metric_table(solver)$name,
               c("leaf_area", "mass_above_ground", "area_stem",
                 "1.net_reproduction_ratio", "1.offspring_production"))
  expect_identical(stand_fitness(solver)[["1.offspring_production"]],
                   solver$patch$offspring_production)
  # The parameters with no column are TF24's, for TF24's reasons.
  expect_identical(census_undifferentiable_tf24f(),
                   census_undifferentiable_tf24())
})

test_that("TF24f's stand gradient is the pinned-step central difference", {
  # The seed (lma), the reproduction schedule (hmat), the stem curve, which
  # moves both ends of the interval and the marginal (stem_P50), and the birth
  # rate. Every column of both metrics agrees to 2e-6; held to 1e-4.
  errors <- tf24f_gradient_errors(
    tf24f_gradient_parameters(), function() Environment("TF24f"),
    c("leaf_area", "1.offspring_production"),
    c("lma", "hmat", "stem_P50", "birth_rate"))
  expect_lt(max(errors), 1e-4)
})

test_that("the gradient holds where the barrier holds the collar at its dry end", {
  # Under dry pulses at k_acclim = 10 the position climbs to the dry end of the
  # interval within a dry spell and the barrier holds it there, while the end
  # itself moves as the soil dries: the collar's rows are then the end's. The
  # stem curve (stem_P50) and the stem's conductance (K_s) move that end; lma
  # does not. Each column agrees to 1e-5; held to 1e-4.
  lifetime <- 10
  errors <- tf24f_gradient_errors(
    tf24f_gradient_parameters(k_acclim = 10, lifetime = lifetime, lma = 0.07),
    function() tf24f_dry_pulses(lifetime),
    "1.offspring_production", c("lma", "stem_P50", "K_s"))
  expect_lt(max(errors), 1e-4)
})

test_that("TF24f's stand gradient approaches TF24's as k_acclim grows", {
  # The forward counterpart is test-model-tf24f.R's convergence test: a large
  # gain recovers TF24's optimum at every step, so the gradient converges too.
  # TF24 is integrated to 1e-6, where its own gradient moves by up to 0.4%
  # against the default tolerance. Measured gap in the gradient of offspring
  # production, the median over the columns below and the largest:
  #
  #   k_acclim      1        10       100
  #   median      0.033    0.0096   0.0010
  #   largest     0.45     0.034    0.0022
  #
  # k = 100 is not run here: the position is stiff at a high gain, so its run
  # takes eight times the steps.
  traits <- paste0("1.", c("lma", "hmat", "vcmax_25", "stem_P50", "K_s",
                           "omega", "k_I", "d_I", "birth_rate"))
  metric <- "1.offspring_production"
  gradient_of <- function(model, k_acclim = 1, ctrl = Control()) {
    solver <- run(tf24f_gradient_parameters(model, k_acclim),
                  Environment(model), ctrl, record = "states")$solver
    g <- stand_gradient(solver, metrics = metric, traits = traits,
                        hyperpar = NULL)
    expect_false(any(stand_gradient_refused(g)))
    g$gradient[metric, ]
  }
  tight <- Control()
  tight$ode_tol_rel <- 1e-6
  tight$ode_tol_abs <- 1e-6
  tf24 <- gradient_of("TF24", ctrl = tight)
  gap <- function(k) abs(gradient_of("TF24f", k) - tf24) / abs(tf24)
  slow <- gap(1)
  fast <- gap(10)
  expect_lt(median(fast), median(slow) / 2)
  expect_lt(max(fast), max(slow))
  expect_lt(median(fast), 0.02)
})
