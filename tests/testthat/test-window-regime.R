## The window regime (#519): a cohort is weighted by when it arrived in one
## long run rather than by the density of patches at its age, so the net
## reproduction ratio is the window average of a cohort's lifetime output.

ff16 <- function(lifetime) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list()
  add_strategies(p, trait_matrix(0.0825, "lma"))
}

patch_of <- function(p) {
  Patch("FF16", "FF16_Env")(p, Environment("FF16"), Control())
}

test_that("a uniform window weights its birth dates equally and nothing outside", {
  p <- window_regime(ff16(30), 10, 20)
  expect_equal(p$patch_type, "window")
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
  p <- window_regime(ff16(30), 10, 20, weight_x = c(10, 20), weight_y = c(1, 3))
  patch <- patch_of(p)
  expect_equal(patch$density(10), 1 / 20)
  expect_equal(patch$density(15), 2 / 20)
  expect_equal(patch$density(20), 3 / 20)
  expect_equal(patch$density(9.99), 0)
  ## Knots inside the window are held constant beyond them.
  p <- window_regime(ff16(30), 10, 20, weight_x = c(12, 14), weight_y = c(2, 2))
  expect_equal(patch_of(p)$density(11), 0.1)
  expect_equal(patch_of(p)$density(19), 0.1)
})

test_that("a window refuses what it cannot weight by", {
  p <- ff16(30)
  expect_error(window_regime(p, 20, 10), "exceed")
  expect_error(window_regime(p, 10, 20, weight_x = c(10, 20), weight_y = 1), "same length")
  expect_error(patch_of(window_regime(p, 10, 20, weight_x = c(10, 20), weight_y = c(1, -1))),
               "non-negative")
  expect_error(patch_of(window_regime(p, 10, 20, weight_x = c(20, 10), weight_y = c(1, 1))),
               "increasing")
  expect_error(patch_of(window_regime(p, 10, 20, weight_x = c(10, 20), weight_y = c(0, 0))),
               "integrate to zero")
  ## The fields default to an empty window, which the type refuses to build.
  q <- ff16(30)
  q$patch_type <- "window"
  expect_error(patch_of(q), "exceed")
})

test_that("the net reproduction ratio is the window average of per-node lifetime output", {
  p <- window_regime(ff16(30), 10, 20)
  scm <- run_scm(p, ctrl = Control(node_density_in_birth_date = TRUE))
  sp <- scm$patch$species[[1]]
  t <- sp$node_times
  fec <- sp$net_reproduction_ratio_by_node
  S_D <- p$strategies[[1]]$pars$S_D
  d <- ifelse(t >= 10 & t <= 20, 1 / 10, 0)
  ## The same trapezia the patch takes, over every node, with the window's
  ## density at each birth date: output times weight over weight alone.
  trapezium <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  expect_equal(scm$net_reproduction_ratios,
               trapezium(t, fec * d * S_D) / trapezium(t, d), tolerance = 1e-12)
  ## Which is the weighted mean of the cohorts' output, so it sits inside their range.
  inside <- t >= 10 & t <= 20
  expect_gte(scm$net_reproduction_ratios, min(fec[inside]) * S_D)
  expect_lte(scm$net_reproduction_ratios, max(fec[inside]) * S_D)
  ## Cohorts born outside the window carry no weight.
  expect_true(all(sp$patch_densities[t < 10 | t > 20] == 0))
  expect_true(all(sp$patch_densities[t >= 10 & t <= 20] == 0.1))
})

test_that("periodic introduction times stop short of the run's end", {
  expect_equal(node_schedule_times_periodic(0, 10, 1), 0:9)
  expect_equal(node_schedule_times_periodic(2.5, 10, 2.5), c(2.5, 5, 7.5))
  expect_error(node_schedule_times_periodic(0, 10, 0), "positive")
  expect_error(node_schedule_times_periodic(10, 0, 1), "exceed")
  p <- ff16(20)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 20, 1))
  scm <- run_scm(p)
  expect_equal(scm$patch$species[[1]]$node_times, 0:19)
})

test_that("the window survives the rebind a gradient runs on", {
  ## rebind_from rebuilds the regime from the parameters it copies; a copy
  ## that left the window behind could not build a "window" patch at all.
  p <- window_regime(ff16(10), 2, 8)
  scm <- run_scm(p, ctrl = Control(node_density_in_birth_date = TRUE),
                 record_trajectory = TRUE)
  g <- stand_gradient(scm, metrics = "leaf_area")
  expect_true(all(is.finite(g$gradient)))
})

test_that("a converged stand's window average does not depend on the window", {
  ## A periodic stand with retirement, run long enough that every cohort in
  ## either window lives out its output before the end. Post hoc, the window
  ## mean of per-node output is what the regime computes (the identity above),
  ## so one run serves every window.
  p <- ff16(500)
  p$patch_type <- "fixed"
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 500, 1))
  scm <- run_scm(p, ctrl = Control(node_density_in_birth_date = TRUE,
                                   retire_density = 1e-8))
  sp <- scm$patch$species[[1]]
  t <- sp$node_times
  out <- sp$net_reproduction_ratio_by_node * p$strategies[[1]]$pars$S_D
  trapezium <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  window_mean <- function(a, b) {
    d <- ifelse(t >= a & t <= b, 1 / (b - a), 0)
    trapezium(t, out * d) / trapezium(t, d)
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
  q <- window_regime(p, 200, 300)
  scm2 <- run_scm(q, ctrl = Control(node_density_in_birth_date = TRUE,
                                    retire_density = 1e-8))
  expect_equal(scm2$net_reproduction_ratios, late, tolerance = 1e-10)
})
