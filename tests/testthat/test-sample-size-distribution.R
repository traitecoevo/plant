## A short two-species FF16 stand with density carried in birth date, which
## sampling needs.
sampling_run <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      p <- scm_base_parameters("FF16")
      p$max_patch_lifetime <- 20
      p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"),
                          birth_rate = c(20, 20))
      cache <<- run_scm(p, ctrl = control(node_density_coordinate = "birth_date"),
                        collect = TRUE)
    }
    cache
  }
})

test_that("sampling refuses a run carried in height, and bad arguments", {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 5
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
  height_run <- run_scm(p, ctrl = control(node_density_coordinate = "height"),
                        collect = TRUE)
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
  times <- c(2, 5, 10, 20)
  area <- 100
  expected <- birth_date_totals(res, times) * area
  set.seed(1)
  x <- sample_size_distribution(res, area = area, times = times)
  observed <- table(factor(x$species, rownames(expected)),
                    factor(x$time, times))
  ## Later times are reached only by following individuals drawn earlier, so
  ## this also checks the thinning by mortality.
  expect_true(all(abs(observed - expected) < 4 * sqrt(expected) + 1))
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
