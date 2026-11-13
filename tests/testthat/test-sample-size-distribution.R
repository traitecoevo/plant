## A short two-species FF16 stand with density carried in birth date, which
## sampling needs.
sampling_run <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      p <- testing_parameters("FF16")
      p$max_patch_lifetime <- 20
      p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"),
                          birth_rate = c(20, 20))
      cache <<- tidy(run(p, control = control(node_density_coordinate = "birth_date"),
                         record = "steps"))
    }
    cache
  }
})

test_that("sampling refuses a run carried in height, and bad arguments", {
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- 5
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
  height_run <- tidy(run(p, control = control(node_density_coordinate = "height"),
                         record = "steps"))
  expect_error(sample_size_distribution(height_run, area = 10), "birth date")

  res <- sampling_run()
  expect_error(sample_size_distribution(res, area = -1), "positive")
  expect_error(sample_size_distribution(res, area = 10, times = 21), "within the run")
  expect_error(sample_size_distribution(res, area = 10, times = c(5, NA)), "finite")
})

test_that("an empty sample has the same columns as a full one", {
  res <- sampling_run()
  set.seed(4)
  full <- sample_size_distribution(res, area = 10)
  empty <- sample_size_distribution(res, area = 1e-12)
  expect_equal(nrow(empty), 0)
  expect_equal(names(empty), names(full))
})

test_that("no individual is born before the first cohort", {
  ## Knots at two requested times for a species first born at t = 5: births in
  ## (1, 5] must not be drawn from a density extrapolated below the first knot.
  knots <- list(
    list(tau = 1, n = 2, mortality = 0, log_density_state = log(2), height = 0.4),
    list(tau = c(5, 7, 10), n = c(2, 2, 2), mortality = c(0.1, 0.1, 0.1),
         log_density_state = log(c(2, 2, 2)), height = c(2, 1, 0.4)))
  set.seed(5)
  x <- sample_species(knots, times = c(1, 10), area = 50, columns = "height")
  expect_gt(nrow(x), 0)
  expect_true(all(x$birth_time >= 5))
})

test_that("the number sampled is Poisson about the integrated density, at every time", {
  res <- sampling_run()
  ## Dense frames, as an animation uses: individuals reach most of them only by
  ## being followed from earlier ones, so this checks the thinning, which must
  ## use the same interpolant the draws do.
  frames <- seq(0.5, 20, by = 0.5)
  check <- c(2, 5, 10, 20)
  ## Large enough that the bias of thinning with a separately interpolated
  ## mortality (about 17% low for species 1 at t = 20) is many SD out.
  area <- 400
  set.seed(1)
  x <- sample_size_distribution(res, area = area, times = frames)
  expected <- birth_date_totals(res, check) * area
  observed <- table(factor(x$species, rownames(expected)),
                    factor(x$time, check))[, as.character(check)]
  expect_true(all(abs(observed - expected) < 4 * sqrt(expected) + 1))

  ## Trees above 2 m, at collected steps, against plant's own integral.
  steps <- vapply(c(10, 20), function(t) {
    res$steps$time[which.min(abs(res$steps$time - t))]
  }, 0)
  set.seed(2)
  y <- sample_size_distribution(res, area = area, times = sort(c(frames, steps)))
  tall <- integrate_over_size_distribution(res$species, lower = 2)
  for (s in rownames(expected)) {
    for (t in steps) {
      e <- tall$density[tall$species == s & tall$time == t] * area
      o <- sum(y$species == s & y$time == t & y$height >= 2)
      expect_lt(abs(o - e), 4 * sqrt(e) + 1)
    }
  }
})

test_that("the birth-date integral agrees with plant's own height integral", {
  res <- sampling_run()
  ## Collected steps, where the height integral is defined.
  times <- vapply(c(5, 10, 20), function(t) {
    res$steps$time[which.min(abs(res$steps$time - t))]
  }, 0)
  by_birth <- birth_date_totals(res, times)
  by_height <- integrate_over_size_distribution(res$species)
  for (s in rownames(by_birth)) {
    for (j in seq_along(times)) {
      h <- by_height$density[by_height$species == s & by_height$time == times[j]]
      expect_equal(by_birth[s, j], h, tolerance = 0.1)
    }
  }
})

test_that("pchip stays within each segment's end values, and integrates exactly", {
  ## Mortality across cohorts at a time when the newest failed to establish:
  ## a monoH.FC spline overshoots far below zero between the first knots.
  x <- c(2.25, 2.5, 3, 3.5, 4)
  y <- c(8.35, 7.62, 6.81, 750, 750)
  f <- pchip(x, y)
  expect_equal(f(x), y)
  for (i in seq_len(length(x) - 1)) {
    v <- f(seq(x[i], x[i + 1], length.out = 101))
    expect_true(all(v >= min(y[i:(i + 1)]) - 1e-9 & v <= max(y[i:(i + 1)]) + 1e-9))
  }
  expect_equal(sum(spline_segments(f, x)),
               stats::integrate(f, min(x), max(x), rel.tol = 1e-10)$value,
               tolerance = 1e-8)
  ## Straight-line data are reproduced exactly.
  g <- pchip(c(0, 1, 3, 4), c(1, 3, 7, 9))
  expect_equal(g(c(0.5, 2, 3.5)), c(2, 5, 8))
})

test_that("individuals persist: consecutive times, growing, never revived", {
  res <- sampling_run()
  times <- seq(1, 20, by = 1)
  set.seed(2)
  x <- sample_size_distribution(res, area = 20, times = times)
  x <- x[order(x$individual, x$time), ]
  same <- diff(x$individual) == 0
  expect_true(all(diff(match(x$time, times))[same] == 1))
  expect_true(all(diff(x$height)[same] >= -1e-9))
  expect_true(all(diff(x$mortality)[same] >= -1e-9))
  expect_true(all(x$birth_time >= 0 & x$birth_time <= x$time))
  ## An individual keeps its species and birth date.
  expect_equal(nrow(unique(x[c("individual", "species", "birth_time")])),
               length(unique(x$individual)))
})

test_that("a sample is reproducible and defaults to a snapshot at the end", {
  res <- sampling_run()
  set.seed(3)
  a <- sample_size_distribution(res, area = 10)
  set.seed(3)
  b <- sample_size_distribution(res, area = 10)
  expect_identical(a, b)
  expect_equal(unique(a$time), max(res$species$time))
  expect_true(all(c("height", "mortality", "area_heartwood") %in% names(a)))
})

test_that("individuals of a cohort that leaves the tracked nodes drop out", {
  ## The oldest node (born at 0) is gone by the second frame, as when a
  ## continuous patch retires a cohort.
  knot <- function(tau) {
    list(tau = tau, n = rep(2, length(tau)), mortality = rep(0.1, length(tau)),
         log_density_state = rep(log(2), length(tau)), height = rev(seq_along(tau)))
  }
  set.seed(6)
  x <- sample_species(list(knot(c(0, 1, 2)), knot(c(1, 2, 3))),
                      times = c(2, 3), area = 50, columns = "height")
  expect_true(any(x$time == 2 & x$birth_time < 1))
  expect_false(any(x$time == 3 & x$birth_time < 1))
})
