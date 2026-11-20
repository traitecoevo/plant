# Fitness as a differentiable stand metric on TF24: each species' net
# reproduction ratio and offspring production, refereed against a central
# difference of runs pinned to the recorded steps, as FF16's are
# (test-ff16-fitness-gradient.R). Each row is refereed on its own scale: the
# fitness sits orders of magnitude away from the census.
#
# The stand is tf24_fitness_stand_parameters(): it reproduces, and is dense
# enough at a birth rate of 10 to hold an equilibrium.

# The stand at a birth rate of 10, recorded once for the two tests that read it.
fitness_stand_at_10 <- local({
  built <- NULL
  function() {
    if (is.null(built)) {
      built <<- tf24_fitness_stand_run(tf24_fitness_stand_parameters(10))
    }
    built
  }
})

test_that("TF24 reports its fitness under the names the gradient uses", {
  solver <- fitness_stand_at_10()
  f <- stand_fitness(solver)
  expect_equal(names(f), c("1.net_reproduction_ratio", "1.offspring_production"))
  expect_identical(f[["1.net_reproduction_ratio"]], solver$patch$net_reproduction_ratios)
  expect_identical(f[["1.offspring_production"]], solver$patch$offspring_production)
  expect_false("S_D" %in% names(census_undifferentiable("TF24")))
  g <- stand_gradient(solver, metrics = names(f), hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  # The birth rate is the last column of the species' block.
  expect_equal(utils::tail(colnames(g$gradient), 1), "1.birth_rate")
  expect_true("1.S_D" %in% colnames(g$gradient))
})

test_that("TF24's fitness and census gradients are the pinned-step central difference", {
  solver <- fitness_stand_at_10()
  fitness <- c("1.net_reproduction_ratio", "1.offspring_production")
  census <- c("leaf_area", "mass_above_ground", "area_stem")
  metrics <- c(fitness, census)
  # Every route into the fecundity state: the reproduction schedule (hmat,
  # a_f1, a_f2), the per-seed cost (omega, a_f3), survival during dispersal,
  # growth and the light field (lma, k_I), mortality (d_I, a_dG1), the leaf
  # (vcmax_25, stem_P50, and the leaf's own R_d_25, jmax_25, root_P50 and
  # TF24_cost_scale, whose rows the leaf supplies), establishment (a_d0) and the
  # birth rate.
  names <- c("hmat", "a_f1", "a_f2", "omega", "a_f3", "S_D", "lma", "k_I",
             "d_I", "a_dG1", "vcmax_25", "stem_P50", "a_d0", "R_d_25",
             "jmax_25", "root_P50", "TF24_cost_scale")
  # The two that reach fecundity alone: the census reads no rate they move, so
  # their census columns are exact zeros rather than differences.
  fecundity_only <- c("a_f3", "S_D")
  g <- stand_gradient(solver, metrics = metrics,
                      traits = paste0("1.", c(names, "birth_rate")),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  expect_gt(g$value[["1.net_reproduction_ratio"]], 1)
  # The branch these rows are taken on: the stand's leaves sit interior.
  kinds <- stats::setNames(census_operating_point_counts(solver, "TF24")[[1]],
                           census_operating_point_names())
  expect_gt(kinds[["interior"]], 0)

  pinned <- pin_steps(solver)
  # lma and omega set the seed height, a root a difference reads at the root
  # tolerance: at offspring_height_tol = 1e-8 their columns fail this test,
  # so it also guards that tolerance. a_d0's difference is flat from 1e-3 to
  # 1e-4 and rounds below it, so it takes 1e-4. Every other column is
  # differenced at 1e-5, where the truncation error on hmat is 6e-5 against
  # 5e-3 at 1e-4. All are held to 1e-4.
  rounding_floor <- "a_d0"
  measure <- function(q) {
    s <- tf24_fitness_stand_run(q, record = FALSE)
    c(stand_fitness(s), stand_census(s))[metrics]
  }
  for (name in c(names, "birth_rate")) {
    rel <- if (name %in% rounding_floor) 1e-4 else 1e-5
    fd <- central_difference(pinned, 1, name, NULL, measure, rel = rel)
    for (metric in metrics) {
      got <- g$gradient[metric, paste0("1.", name)]
      if (metric %in% census && name %in% fecundity_only) {
        expect_identical(got, 0, label = paste(metric, name))
      } else {
        expect_lt(abs(got - fd[[metric]]), 1e-4 * abs(fd[[metric]]),
                  label = paste(metric, name))
      }
    }
  }
  # S_D multiplies the whole integral, so its column is the value over S_D.
  S_D <- pinned$strategies[[1]]$pars$S_D
  expect_equal(g$gradient["1.net_reproduction_ratio", "1.S_D"],
               g$value[["1.net_reproduction_ratio"]] / S_D, tolerance = 1e-12)
})

test_that("Newton on d ln J / d b finds the equilibrium birth rate, and an invader is differentiated there", {
  # The resident at its demographic equilibrium makes as many seeds as it is
  # given: J(b*) = b*. Newton in ln b on ln J - ln b, with the slope from the
  # birth-rate column of offspring production J.
  b <- 10
  for (k in 1:8) {
    resident <- tf24_fitness_stand_run(tf24_fitness_stand_parameters(b))
    g <- stand_gradient(resident, metrics = "1.offspring_production",
                        traits = "1.birth_rate", hyperpar = NULL)
    J <- g$value[["1.offspring_production"]]
    excess <- log(J) - log(b)
    if (abs(excess) < 1e-8) {
      break
    }
    b <- exp(log(b) - excess / (b * g$gradient[[1]] / J - 1))
  }
  # Linear, not quadratic: the adaptive runs' slope in b differs from the
  # swept slope of their own steps by about 2.5%. Seven runs here.
  expect_lt(abs(excess), 1e-8)
  expect_equal(resident$patch$net_reproduction_ratios, 1, tolerance = 1e-8)

  # An invader at the resident's traits has the resident's ratio, one.
  resident$keep_states <- TRUE
  resident$run_mutant(tf24_fitness_stand_parameters(1))
  expect_equal(resident$patch$net_reproduction_ratios, 1, tolerance = 1e-8)

  # An invader at lma x0.9, differentiated in the recorded field of the
  # resident's equilibrium: the sweep against a central difference of
  # invasions of the same recording.
  invader <- tf24_fitness_stand_parameters(1, 0.9)
  resident$run_mutant(invader)
  gi <- stand_gradient(resident, metrics = "1.net_reproduction_ratio",
                       traits = c("1.lma", "1.hmat"), hyperpar = NULL)
  expect_false(any(stand_gradient_refused(gi)))
  expect_gt(gi$value[["1.net_reproduction_ratio"]], 1)
  resident$keep_states <- FALSE
  for (name in c("lma", "hmat")) {
    x <- invader$strategies[[1]]$pars[[name]]
    h <- 1e-4 * x
    at <- function(v) {
      q <- invader
      q$strategies[[1]]$pars[[name]] <- v
      resident$run_mutant(q)
      resident$patch$net_reproduction_ratios
    }
    fd <- (at(x + h) - at(x - h)) / (2 * h)
    expect_equal(gi$gradient["1.net_reproduction_ratio", paste0("1.", name)],
                 fd, tolerance = if (name == "lma") 1e-3 else 1e-4, label = name)
  }
})
