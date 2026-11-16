# TF24f's stand gradient: TF24's census and fitness, with the tracked collar
# potential carried through the sweep as an ODE state. Refereed against central
# differences of runs pinned to the recorded steps, as TF24's and FF16's are
# (test-tf24-fitness-gradient.R, test-ff16-stand-gradient.R).
#
# The stand is the one test-model-tf24f.R runs: five years, hmat 5, a birth rate
# of 20, on the birth-date coordinate. There the tracked collar stays inside its
# interval at every rate evaluation and the leaf never shuts, so this referees
# the tracked state's own rows: the leaf billed at the state, the marginal's
# trait, soil and light rows, its curvature, and the birth seed.

tf24f_gradient_parameters <- function(model = "TF24f", k_acclim = 1) {
  p <- testing_parameters(model)
  p$max_patch_lifetime <- 5
  p <- add_strategies(p, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                      hyperpar = default_hyperpar(model), birth_rate = 20)
  if (model == "TF24f") {
    s <- p$strategies[[1]]
    s$k_acclim <- k_acclim
    p$strategies[[1]] <- s
  }
  p
}

tf24f_gradient_run <- function(p, ctrl = Control(), record = TRUE) {
  run(p, Environment("TF24f"), ctrl,
      record = if (record) "states" else character(0))$solver
}

tf24f_gradient_metrics <- c("leaf_area", "1.offspring_production")

tf24f_gradient_of <- function(solver, traits) {
  stand_gradient(solver, metrics = tf24f_gradient_metrics,
                 traits = paste0("1.", traits), hyperpar = NULL)
}

test_that("TF24f declares TF24's census and fitness", {
  solver <- tf24f_gradient_run(tf24f_gradient_parameters())
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
  solver <- tf24f_gradient_run(tf24f_gradient_parameters())
  # The seed (lma), the reproduction schedule (hmat), the leaf through the
  # marginal (vcmax_25) and the stem curve (stem_P50, K_s), and the birth rate.
  # At a relative step of 1e-5 every column agrees to 4e-6.
  traits <- c("lma", "hmat", "vcmax_25", "stem_P50", "K_s", "birth_rate")
  g <- tf24f_gradient_of(solver, traits)
  expect_false(any(stand_gradient_refused(g)))

  pinned <- solver$parameters
  pinned$ode_times <- solver$ode_times
  pinned$ode_step_sizes <- solver$ode_step_sizes
  metrics_of <- function(q) {
    s <- tf24f_gradient_run(q, record = FALSE)
    c(stand_census(s)["leaf_area"], stand_fitness(s)["1.offspring_production"])
  }
  expect_equal(metrics_of(pinned), g$value[tf24f_gradient_metrics],
               tolerance = 1e-12)
  for (name in traits) {
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
    h <- 1e-5 * abs(x)
    fd <- (metrics_of(set(x + h)) - metrics_of(set(x - h))) / (2 * h)
    for (m in tf24f_gradient_metrics) {
      expect_lt(abs(g$gradient[m, paste0("1.", name)] - fd[[m]]),
                1e-4 * abs(fd[[m]]), label = paste(m, name))
    }
  }
})

test_that("TF24f's stand gradient approaches TF24's as k_acclim grows", {
  # The forward counterpart is test-model-tf24f.R's convergence test: a large
  # gain recovers TF24's optimum at every step, so the gradient converges too.
  # TF24 is integrated to 1e-6, where its own gradient moves by up to 0.4%
  # against the default tolerance. Measured gap in offspring production, the
  # median over the columns below and the largest:
  #
  #   k_acclim      1       10      100
  #   median      0.213   0.058   0.0065
  #   largest     4.30    0.545   0.066
  #
  # k = 100 is not run here: its sweep costs about twelve seconds, against
  # three at k = 1, because the tracked state is stiff at a high gain.
  traits <- c("lma", "hmat", "vcmax_25", "stem_P50", "K_s", "omega", "k_I",
              "d_I", "birth_rate")
  tight <- Control()
  tight$ode_tol_rel <- 1e-6
  tight$ode_tol_abs <- 1e-6
  tf24 <- tf24f_gradient_of(
    tf24f_gradient_run(tf24f_gradient_parameters("TF24"), tight), traits)
  gap <- function(k) {
    g <- tf24f_gradient_of(
      tf24f_gradient_run(tf24f_gradient_parameters(k_acclim = k)), traits)
    expect_false(any(stand_gradient_refused(g)))
    abs(g$gradient["1.offspring_production", ] -
          tf24$gradient["1.offspring_production", ]) /
      abs(tf24$gradient["1.offspring_production", ])
  }
  slow <- gap(1)
  fast <- gap(10)
  expect_lt(median(fast), median(slow) / 2)
  expect_lt(max(fast), max(slow))
  expect_lt(median(fast), 0.1)
})
