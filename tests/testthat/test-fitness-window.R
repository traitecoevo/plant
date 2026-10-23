## The fitness window (#519): a cohort is weighted by when it arrived in one
## long run rather than by the density of patches at its age, so the net
## reproduction ratio is the window average of a cohort's lifetime output.

ff16 <- function(lifetime) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list()
  add_strategies(p, trait_matrix(0.0825, "lma"))
}

## The regime's ratio, replicated: over the span of cohorts the window weights,
## output times weight by establishment weight, over the trapezium of the weight.
window_mean_of <- function(sp, out, d) {
  span <- range(which(d > 0))
  k <- span[1]:span[2]
  t <- sp$node_times[k]
  trapezium <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  sum(establishment_weights_within(sp, k) * out[k] * d[k]) / trapezium(t, d[k])
}

patch_of <- function(p) {
  Patch("FF16", "FF16_Env")(p, Environment("FF16"), Control())
}

test_that("a uniform window weights its birth dates equally and nothing outside", {
  p <- fitness_window(ff16(30), 10, 20)
  expect_equal(p$patch_type, "fitness_window")
  patch <- patch_of(p)
  expect_equal(patch$density(5), 0)
  expect_equal(patch$density(10), 0.1)
  expect_equal(patch$density(15), 0.1)
  expect_equal(patch$density(20), 0.1)
  expect_equal(patch$density(25), 0)
  expect_equal(patch$pr_survival(15), 1)
  expect_equal(patch$pr_survival(25), 1)
})

test_that("a tabulated weight is read at the birth date and normalised over the window", {
  ## Weight 1 at 10 rising to 3 at 20: integral 20, so the density runs from
  ## 1/20 to 3/20 and is 2/20 in the middle.
  p <- fitness_window(ff16(30), 10, 20, weight_x = c(10, 20), weight_y = c(1, 3))
  patch <- patch_of(p)
  expect_equal(patch$density(10), 1 / 20)
  expect_equal(patch$density(15), 2 / 20)
  expect_equal(patch$density(20), 3 / 20)
  expect_equal(patch$density(9.99), 0)
  ## Knots inside the window are held constant beyond them.
  p <- fitness_window(ff16(30), 10, 20, weight_x = c(12, 14), weight_y = c(2, 2))
  expect_equal(patch_of(p)$density(11), 0.1)
  expect_equal(patch_of(p)$density(19), 0.1)
})

test_that("a window refuses what it cannot weight by", {
  p <- ff16(30)
  expect_error(fitness_window(p, 20, 10), "exceed")
  expect_error(fitness_window(p, 10, 20, weight_x = c(10, 20), weight_y = 1), "same length")
  expect_error(patch_of(fitness_window(p, 10, 20, weight_x = c(10, 20), weight_y = c(1, -1))),
               "non-negative")
  expect_error(patch_of(fitness_window(p, 10, 20, weight_x = c(20, 10), weight_y = c(1, 1))),
               "increasing")
  expect_error(patch_of(fitness_window(p, 10, 20, weight_x = c(10, 20), weight_y = c(0, 0))),
               "integrate to zero")
  ## The fields default to an empty window, which the type refuses to build.
  q <- ff16(30)
  q$patch_type <- "fitness_window"
  expect_error(patch_of(q), "exceed")
})

test_that("the net reproduction ratio is the window average of per-node lifetime output", {
  p <- fitness_window(ff16(30), 10, 20)
  scm <- run_scm(p, ctrl = Control(node_density_in_birth_date = TRUE))
  sp <- scm$patch$species[[1]]
  t <- sp$node_times
  fec <- sp$net_reproduction_ratio_by_node
  S_D <- p$strategies[[1]]$pars$S_D
  d <- ifelse(t >= 10 & t <= 20, 1 / 10, 0)
  ## The same sum the patch takes, over the window's cohorts, with the
  ## window's density at each birth date: output times weight over weight.
  expect_equal(scm$net_reproduction_ratios, window_mean_of(sp, fec * S_D, d),
               tolerance = 1e-12)
  ## Which is the mean output per arriving seed, so no more than the most a
  ## seed that establishes produces.
  inside <- t >= 10 & t <= 20
  expect_gt(scm$net_reproduction_ratios, 0)
  expect_lte(scm$net_reproduction_ratios, max(fec[inside]) * S_D)
  ## Cohorts born outside the window carry no weight.
  expect_true(all(sp$patch_densities[t < 10 | t > 20] == 0))
  expect_true(all(sp$patch_densities[t >= 10 & t <= 20] == 0.1))
})

## A guard for any change to the fitness integral (#658 routes it through a
## plain trapezium): offspring production takes the same normalised window mean,
## of output scaled by the birth rate at each birth date.
test_that("offspring production is the window mean of output times birth rate", {
  p <- fitness_window(ff16(30), 10, 20)
  p$strategies[[1]]$birth_rate_y <- 3
  scm <- run_scm(p, ctrl = Control(node_density_in_birth_date = TRUE))
  expect_equal(scm$offspring_production, 3 * scm$net_reproduction_ratios,
               tolerance = 1e-12)
})

test_that("a window the run cannot score is refused before the run", {
  ## Past the run's end, its last arrivals are never integrated.
  expect_error(run_scm(fitness_window(ff16(30), 10, 30)), "run's end")
  ## With no cohort inside it, there is nothing to score.
  p <- fitness_window(ff16(30), 10, 20)
  p$node_schedule_times <- list(c(0, 1, 2, 25))
  expect_error(run_scm(p), "No cohort")
  ## And a misspelt patch type is refused rather than run unweighted.
  p <- fitness_window(ff16(30), 10, 20)
  p$patch_type <- "fitness-window"
  expect_error(run_scm(p), "Unknown patch_type")
})

test_that("periodic introduction times stop short of the run's end", {
  expect_equal(node_schedule_times_periodic(0, 10, 1), 0:9)
  expect_equal(node_schedule_times_periodic(2.5, 10, 2.5), c(2.5, 5, 7.5))
  # A window weighs a cohort at its end, so an invader's schedule can keep it.
  expect_equal(node_schedule_times_periodic(2.5, 10, 2.5, include_end = TRUE),
               c(2.5, 5, 7.5, 10))
  expect_error(node_schedule_times_periodic(0, 10, 0), "positive")
  expect_error(node_schedule_times_periodic(10, 0, 1), "exceed")
  p <- ff16(20)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 20, 1))
  scm <- run_scm(p)
  expect_equal(scm$patch$species[[1]]$node_times, 0:19)
})

test_that("the window survives the rebind a gradient runs on", {
  ## rebind_from rebuilds the regime from the parameters it copies; a copy
  ## that left the window behind could not build a "fitness_window" patch at all.
  p <- fitness_window(ff16(10), 2, 8)
  scm <- run_scm(p, ctrl = Control(node_density_in_birth_date = TRUE),
                 record_trajectory = TRUE)
  g <- stand_gradient(scm, metrics = "leaf_area", hyperpar = NULL)
  expect_true(all(is.finite(g$gradient)))
})

test_that("a converged stand's window average does not depend on the window", {
  ## A periodic stand run long enough that every cohort in either window lives
  ## out its output before the end. Post hoc, the window mean of per-node output
  ## is what the regime computes (the identity above), so one run serves every
  ## window. No top layer here: it would absorb the window's cohorts in the run
  ## without a window, stopping their per-node output where it took them (the
  ## rest goes to the layer), and the identity is per node.
  p <- ff16(500)
  p$patch_type <- "fixed"
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 500, 1))
  scm <- run_scm(p, ctrl = Control(node_density_in_birth_date = TRUE))
  sp <- scm$patch$species[[1]]
  t <- sp$node_times
  out <- sp$net_reproduction_ratio_by_node * p$strategies[[1]]$pars$S_D
  window_mean <- function(a, b) {
    window_mean_of(sp, out, ifelse(t >= a & t <= b, 1 / (b - a), 0))
  }
  early <- window_mean(100, 200)
  late <- window_mean(200, 300)
  long <- window_mean(100, 300)
  expect_equal(late, early, tolerance = 1e-3)
  expect_equal(long, early, tolerance = 1e-3)
  ## The founding century is not converged, and a window that ends too close
  ## to the run's end truncates its last arrivals: both read differently.
  expect_gt(abs(window_mean(0, 100) / early - 1), 1e-2)
  expect_gt(abs(window_mean(400, 500) / early - 1), 1e-2)
  ## And the regime run with one of these windows agrees with the post-hoc mean.
  q <- fitness_window(p, 200, 300)
  scm2 <- run_scm(q, ctrl = Control(node_density_in_birth_date = TRUE))
  expect_equal(scm2$net_reproduction_ratios, late, tolerance = 1e-10)
})

test_that("cohorts seeded from a burn-in weigh nothing before the window", {
  ## A burn-in under no regime exports its stand; resumed into a window, the
  ## seeded cohorts were born before it and carry no weight, whatever the
  ## burn-in stamped them with, so the ratio is the mean over the window's own
  ## arrivals.
  p <- ff16(20)
  p$patch_type <- "fixed"
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 20, 1))
  ctrl <- Control(node_density_in_birth_date = TRUE)
  burn <- run_scm(p, ctrl = ctrl)
  state <- export_patch_state(burn)
  expect_true(all(state$patch_density[[1]] == 1))
  p2 <- fitness_window(ff16(40), 25, 35)
  p2 <- set_initial_state(p2, state)
  p2$node_schedule_times <- list(node_schedule_times_periodic(20, 40, 1))
  scm <- run_scm(p2, ctrl = ctrl)
  sp <- scm$patch$species[[1]]
  t <- sp$node_times
  expect_true(all(sp$patch_densities[t < 25] == 0))
  expect_true(all(sp$patch_densities[t >= 25 & t <= 35] == 0.1))
  out <- sp$net_reproduction_ratio_by_node * p2$strategies[[1]]$pars$S_D
  d <- ifelse(t >= 25 & t <= 35, 0.1, 0)
  expect_equal(scm$net_reproduction_ratios, window_mean_of(sp, out, d),
               tolerance = 1e-12)
  ## A resume under the regime that stamped the export is unchanged by this.
  p3 <- set_initial_state(ff16(40), state)
  p3$patch_type <- "fixed"
  p3$node_schedule_times <- list(node_schedule_times_periodic(20, 40, 1))
  scm3 <- run_scm(p3, ctrl = ctrl)
  expect_true(all(scm3$patch$species[[1]]$patch_densities == 1))
})

test_that("cohorts outside the window do not move the ratio", {
  ## A resident whose schedule starts before the window, and a copy of the same
  ## strategy a million times rarer whose schedule starts at the window, in one
  ## run. Both stand in the same field with the same strategy and birth dates,
  ## and a cohort's per-capita output does not depend on its own density, so
  ## their window cohorts' outputs are the same numbers and the two must score
  ## the same. Over every node, the resident's zero-weight cohort before the
  ## window would give its first window cohort a whole interval and the copy's
  ## half of one. Long enough for the window's cohorts to reproduce, or both
  ## ratios are near zero and agree vacuously.
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 100
  p$node_schedule_times <- list()
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.0825), "lma"),
                      birth_rate = list(1, 1e-6))
  p <- fitness_window(p, 10, 20)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 100, 1),
                                node_schedule_times_periodic(10, 100, 1))
  scm <- run_scm(p, ctrl = Control(node_density_in_birth_date = TRUE))
  r <- scm$net_reproduction_ratios
  expect_gt(r[[1]], 1)
  expect_equal(r[[2]], r[[1]], tolerance = 1e-12)
})

test_that("a tabulated weight in a full run is the weighted mean of the cohorts", {
  ## A weight that rises and falls over the window, read at each cohort's birth
  ## date: the run's ratio is the weighted mean of the per-node outputs of the
  ## same run without the window, which the window does not change.
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- 120
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 120, 1))
  ctrl <- Control(node_density_in_birth_date = TRUE)
  plain <- run_scm(p, ctrl = ctrl)
  sp <- plain$patch$species[[1]]
  t <- sp$node_times
  out <- sp$net_reproduction_ratio_by_node * p$strategies[[1]]$pars$S_D
  x <- c(20, 30, 40)
  y <- c(1, 3, 1)
  w <- ifelse(t >= 20 & t <= 40, approx(x, y, t, rule = 2)$y, 0)
  windowed <- run_scm(fitness_window(p, 20, 40, weight_x = x, weight_y = y), ctrl = ctrl)
  expect_gt(windowed$net_reproduction_ratios, 1)
  expect_equal(windowed$net_reproduction_ratios, window_mean_of(sp, out, w),
               tolerance = 1e-10)
  ## And it is not the uniform mean.
  expect_false(isTRUE(all.equal(windowed$net_reproduction_ratios,
                                window_mean_of(sp, out, as.numeric(w > 0)),
                                tolerance = 1e-6)))
})
