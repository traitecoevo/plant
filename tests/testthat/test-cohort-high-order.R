# control$cohort_high_order: product integration for the competition field and
# Simpson's rule for lifetime fitness, in place of the trapezium.

test_that("Simpson's rule is exact for quadratics on a non-uniform grid", {
  set.seed(1)
  f <- function(x) 3 - 2 * x + 5 * x^2
  F <- function(x) 3 * x - x^2 + 5 * x^3 / 3
  for (n in c(3, 4, 7, 10)) {
    x <- sort(runif(n, 0, 2))
    expect_equal(simpson(x, f(x)), F(max(x)) - F(min(x)), tolerance = 1e-13)
  }
  # Two points take the trapezium.
  expect_equal(simpson(c(0, 2), c(1, 3)), 4)
  # A zero-width interval leaves a pair without a parabola, and the pair takes
  # the trapezium on its two intervals.
  expect_equal(simpson(c(0, 0, 1), c(5, 1, 3)), 2)
})

# One FF16 species, with the rules on or off.
high_order_run <- function(on, schedule = NULL, refine = FALSE) {
  p <- testing_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  if (!is.null(schedule)) {
    p$node_schedule_times <- list(schedule)
  }
  ctrl <- control()
  ctrl$cohort_high_order <- on
  run(p, control = ctrl, refine_schedule = refine)
}

test_that("The field build agrees with the per-height reduction", {
  r <- high_order_run(TRUE)
  pa <- r$solver$patch
  st <- pa$environment$light_availability$state
  ok <- st[, "light_availability"] < 1
  field <- -log(st[ok, "light_availability"])
  walk <- vapply(st[ok, "height"], function(z) pa$compute_competition(z), 0.0)
  # The trapezium's own field and walk agree to this level too: it is reading
  # the field back through the stored light availability.
  expect_equal(field, walk, tolerance = 1e-6)
})

test_that("Fitness is Simpson's rule over the cohorts' birth dates", {
  r <- high_order_run(TRUE)
  sp <- r$solver$patch$species[[1]]
  y <- sp$net_reproduction_ratio_by_node * sp$patch_densities *
    r$scenario$p$strategies[[1]]$pars$S_D
  expect_equal(r$net_reproduction_ratios[[1]], simpson(sp$node_times, y),
               tolerance = 1e-13)
  expect_false(isTRUE(all.equal(r$net_reproduction_ratios[[1]],
                                trapezium(sp$node_times, y), tolerance = 1e-6)))
})

test_that("The high-order rules are closer to a fine schedule's answer", {
  # The schedule a trapezium run refines to, and the reference on it bisected
  # three times, which agrees with twice bisected to 3e-5.
  coarse <- high_order_run(FALSE, refine = TRUE)$scenario$p$node_schedule_times[[1]]
  fine <- coarse
  for (k in 1:3) {
    fine <- sort(c(fine, head(fine, -1) + diff(fine) / 2))
  }
  heights <- seq(0, 15, by = 1)
  light <- function(r) {
    vapply(heights, function(z) r$solver$patch$compute_competition(z), 0.0)
  }
  ref <- high_order_run(TRUE, fine)
  err <- function(on) {
    r <- high_order_run(on, coarse)
    c(fitness = abs(r$net_reproduction_ratios[[1]] /
                      ref$net_reproduction_ratios[[1]] - 1),
      light = max(abs(light(r) / light(ref) - 1)))
  }
  trap <- err(FALSE)
  high <- err(TRUE)
  # On this stand the high-order errors are about a quarter (fitness) and a
  # sixth (light) of the trapezium's; the test asks for half.
  expect_lt(high[["fitness"]], trap[["fitness"]] / 2)
  expect_lt(high[["light"]], trap[["light"]] / 2)
})

test_that("The high-order rules refuse what they do not cover", {
  ctrl <- control()
  ctrl$cohort_high_order <- TRUE
  # TF24 carries its density in birth date.
  p <- testing_parameters("TF24")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  expect_error(run(p, control = ctrl), "height size-density coordinate")
  # The absorbing top layer.
  p <- testing_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  ctrl_top <- ctrl
  ctrl_top$node_density_coordinate <- "height"
  ctrl_top$absorb_epsilon <- 0.01
  expect_error(run(p, control = ctrl_top), "absorbing top layer")
  # The census, whose reductions take the trapezium.
  r <- run(p, control = ctrl)
  expect_error(stand_census(r), "control\\$cohort_high_order is on")
})
