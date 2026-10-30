# Pulse-series generators (#522, #519): a series is built up front, as data, and
# handed to a model's pulse constructor.

test_that("a drawn series is a stationary marked Poisson process", {
  set.seed(1)
  s <- pulse_series_poisson(0, 1000, rate = 20, mean_amount = 0.014)
  expect_named(s, c("time", "amount"))
  expect_false(is.unsorted(s$time))
  expect_true(all(s$time > 0 & s$time <= 1000))
  expect_true(all(s$amount > 0))
  # 20,000 pulses expected, sd about 141; the mean amount to about 1 %.
  expect_lt(abs(nrow(s) - 20000), 4 * sqrt(20000))
  expect_equal(mean(s$amount), 0.014, tolerance = 0.02)
})

test_that("a drawn series draws count, then times, then amounts", {
  # The order is the contract that lets a script's series be rebuilt from its
  # seed, and the #519 pilot's committed numbers rest on it.
  set.seed(7)
  s <- pulse_series_poisson(0, 50, rate = 3, mean_amount = 0.02)
  set.seed(7)
  n <- rpois(1, 3 * 50)
  time <- sort(runif(n, 0, 50))
  amount <- rexp(n, rate = 1 / 0.02)
  expect_identical(s$time, time)
  expect_identical(s$amount, amount)
})

test_that("a regular series puts one pulse every interval", {
  s <- pulse_series_regular(2, 3, interval = 0.25, amount = 0.01)
  expect_equal(s$time, c(2.25, 2.5, 2.75, 3))
  expect_equal(s$amount, rep(0.01, 4))
  ## A last pulse landing on `to` is on it, not a rounding error past it, so a
  ## run whose horizon is `to` accepts it.
  r <- pulse_series_regular(0, 0.7, interval = 0.1, amount = 0.01)
  expect_true(all(r$time <= 0.7))
  expect_identical(r$time[[nrow(r)]], 0.7)
})

test_that("a record repeats end to end over a longer run", {
  s <- pulse_series_record(time = c(0.1, 0.5, 0.9), amount = c(1, 2, 3),
                           from = 10, to = 12.5)
  # The span is (from, to], so the pulse landing on `to` itself is in it.
  expect_equal(s$time, c(10.1, 10.5, 10.9, 11.1, 11.5, 11.9, 12.1, 12.5))
  expect_equal(s$amount, c(1, 2, 3, 1, 2, 3, 1, 2))
  # Over one period it is the record itself, shifted to start at `from`.
  one <- pulse_series_record(time = c(0.1, 0.5, 0.9), amount = c(1, 2, 3))
  expect_equal(one$time, c(0.1, 0.5, 0.9))
})

test_that("generators refuse what is not a series", {
  expect_error(pulse_series_poisson(5, 5, 1, 1), "to")
  expect_error(pulse_series_poisson(0, 5, -1, 1), "rate")
  expect_error(pulse_series_regular(0, 5, 0, 1), "interval")
  expect_error(pulse_series_record(c(0.5, 0.2), c(1, 1)), "increasing")
  expect_error(pulse_series_record(c(0, 0.5), c(1, 1)), "\\(0, period\\]")
  expect_error(pulse_series_record(c(0.5, 2), c(1, 1), period = 1), "\\(0, period\\]")
})

test_that("a series becomes rainfall pulses", {
  set.seed(2)
  s <- pulse_series_poisson(0, 2, rate = 10, mean_amount = 0.01)
  p <- add_strategies(scm_base_parameters("TF24"), trait_matrix(0.0825, "lma"),
                      hyperpar = TF24_hyperpar)
  ev <- events(events_default(p), rainfall_pulse(s$time, s$amount))
  pulses <- ev$type == "resource_pulse"
  expect_equal(sum(pulses), nrow(s))
  expect_equal(ev$time[pulses], s$time)
})

test_that("a series drawn from a fitted record keeps its total and its season", {
  set.seed(3)
  date <- seq(as.Date("1981-01-01"), as.Date("2020-12-31"), by = "day")
  summer <- format(date, "%m") %in% c("12", "01", "02")
  amount <- ifelse(runif(length(date)) < ifelse(summer, 0.5, 0.05),
                   rexp(length(date), rate = 1 / 0.01), 0)
  fit <- pulse_series_fit(date, amount)
  expect_equal(fit$years, 1981:2020)

  set.seed(4)
  s <- pulse_series_markov(0, 200, fit)
  expect_false(is.unsorted(s$time, strictly = TRUE))
  expect_true(all(s$time > 0 & s$time <= 200 & s$amount > 0))
  # The mean annual total of 200 drawn years, about 0.6 m with a CV near 0.2,
  # is the record's to a few percent.
  annual <- tapply(amount, format(date, "%Y"), sum)
  expect_equal(sum(s$amount) / 200, mean(annual), tolerance = 0.05)
  # January has ten times July's chance of rain.
  day <- floor((s$time %% 1) * 365)
  expect_gt(sum(day < 31), 5 * sum(day >= 181 & day < 212))

  # Pooling days into weekly pulses moves the water, not its amount.
  set.seed(4)
  w <- pulse_series_markov(0, 200, fit, days_per_pulse = 7)
  expect_equal(sum(w$amount), sum(s$amount))
  expect_true(all(w$time > 0 & w$time <= 200))
})

test_that("a fit refuses a record with no complete year", {
  date <- seq(as.Date("2001-03-01"), as.Date("2002-02-28"), by = "day")
  expect_error(pulse_series_fit(date, rep(0.001, length(date))), "complete year")
})
