## Pulse-series generators (#522, #519). The decision on #522 is that the whole
## series is baked in up front -- drawn, regular, or from an observed record --
## and handed to the run as events, so a series is plain data the run, a mutant
## and a gradient all see the same of. These build the series; a model's own
## constructor turns it into events (TF24: rainfall_pulse(s$time, s$amount)).
## Nothing here is specific to a model or to water.

##' Series of resource pulses
##'
##' Three ways to build a series of pulses for a run, each returning the pulses'
##' times and amounts as a data frame, sorted by time, to hand to a model's pulse
##' constructor: for TF24, \code{rainfall_pulse(s$time, s$amount)}.
##'
##' \code{pulse_series_poisson()} draws a stationary marked Poisson process:
##' pulses arrive at \code{rate} per year, and each one's amount is exponential
##' with mean \code{mean_amount}, so the long-run supply is \code{rate *
##' mean_amount} per year. It draws from R's random number stream, so
##' \code{set.seed()} fixes the series, and one realisation shared by every run
##' that is compared is what makes a comparison of fitness a comparison of
##' traits (common random numbers). It draws the count, then the times, then the
##' amounts, so a series over a longer span does not contain a shorter one.
##'
##' \code{pulse_series_regular()} puts a pulse of \code{amount} every
##' \code{interval} years.
##'
##' \code{pulse_series_record()} takes an observed record -- each wet day a
##' pulse, say, with times in years -- and, given \code{to} beyond the record's
##' end, repeats it end to end to cover a longer run: a 100-year record behind a
##' burn-in and a window of several centuries. Repeating a record makes the
##' series periodic, which a window shorter than the record does not see.
##'
##' @title Series of resource pulses
##' @param from,to The span the series covers: pulses fall in \code{(from, to]}.
##' @param rate Mean number of pulses per year.
##' @param mean_amount Mean amount per pulse, in the units the pulse constructor
##'   takes (TF24: metres of water).
##' @return A data frame with columns \code{time} and \code{amount}.
##' @seealso \code{\link{rainfall_pulse}}, \code{\link{events}}
##' @export
##' @examples
##' set.seed(1)
##' s <- pulse_series_poisson(0, 100, rate = 20, mean_amount = 0.014)
##' nrow(s)                      # about 2000 pulses
##' sum(s$amount) / 100          # about 0.28 m of water a year
##' p <- add_strategies(scm_base_parameters("TF24"), trait_matrix(0.0825, "lma"))
##' ev <- events(events_default(p), rainfall_pulse(s$time, s$amount))
pulse_series_poisson <- function(from, to, rate, mean_amount) {
  check_pulse_span(from, to)
  check_positive_scalar(rate, "rate")
  check_positive_scalar(mean_amount, "mean_amount")
  n <- stats::rpois(1, rate * (to - from))
  time <- sort(stats::runif(n, from, to))
  amount <- stats::rexp(n, rate = 1 / mean_amount)
  ## runif can return its lower bound in principle; (from, to] is the contract.
  keep <- time > from
  data.frame(time = time[keep], amount = amount[keep])
}

##' @param interval Years between pulses.
##' @param amount Amount per pulse (\code{pulse_series_regular}), or per
##'   recorded pulse (\code{pulse_series_record}).
##' @rdname pulse_series_poisson
##' @export
pulse_series_regular <- function(from, to, interval, amount) {
  check_pulse_span(from, to)
  check_positive_scalar(interval, "interval")
  check_positive_scalar(amount, "amount")
  n <- floor((to - from) / interval + 1e-10)
  # The slack that keeps a pulse landing on `to` can put it a rounding error
  # past it, outside (from, to]; it belongs on `to`.
  time <- pmin(from + seq_len(n) * interval, to)
  data.frame(time = time, amount = rep(amount, length(time)))
}

##' @param time Times of the recorded pulses, in years from the start of the
##'   record, in \code{(0, period]} and increasing: day of record / 365.25, say.
##' @param period The record's length in years, the spacing at which it
##'   repeats. Defaults to its last time rounded up to a whole year.
##' @rdname pulse_series_poisson
##' @export
pulse_series_record <- function(time, amount, from = 0, to = from + period,
                                period = ceiling(max(time))) {
  if (length(time) == 0L || length(time) != length(amount)) {
    stop("`time` and `amount` must be non-empty and the same length")
  }
  if (any(!is.finite(time)) || is.unsorted(time, strictly = TRUE)) {
    stop("`time` must be finite and strictly increasing")
  }
  if (any(!is.finite(amount)) || any(amount < 0)) {
    stop("`amount` must be finite and non-negative")
  }
  check_positive_scalar(period, "period")
  if (time[[1]] <= 0 || time[[length(time)]] > period) {
    stop("`time` must lie in (0, period]: times run from the start of the record")
  }
  check_pulse_span(from, to)
  n_repeat <- ceiling((to - from) / period)
  t_all <- from + as.vector(outer(time, (seq_len(n_repeat) - 1) * period, `+`))
  a_all <- rep(amount, n_repeat)
  keep <- t_all <= to
  data.frame(time = t_all[keep], amount = a_all[keep])
}

check_pulse_span <- function(from, to) {
  if (!(length(from) == 1L && length(to) == 1L && is.finite(from) &&
        is.finite(to) && to > from)) {
    stop("`from` and `to` must be finite scalars with `to` > `from`")
  }
}

check_positive_scalar <- function(x, name) {
  if (!(length(x) == 1L && is.finite(x) && x > 0)) {
    stop("`", name, "` must be a positive finite scalar")
  }
}

##' Series of resource pulses fitted to a daily record
##'
##' \code{pulse_series_fit()} fits a seasonal weather generator to an observed
##' daily record, and \code{pulse_series_markov()} draws a series of any length
##' from the fit. Where \code{pulse_series_record()} repeats a record, which
##' makes the series periodic, a fitted series is stationary, and a different
##' realisation for each seed.
##'
##' The fit uses the record's complete years only: those with an observation on
##' every day. Each year is put in a class by its total -- dry, normal or wet,
##' split at the terciles of the annual totals. For each month and class, the
##' fit keeps the chance that a day is wet after a wet day and after a dry one
##' (a first-order Markov chain), and the empirical quantiles of the amounts on
##' wet days. A month and class with fewer than \code{min_wet} wet days, or
##' fewer than \code{min_wet} days after a wet or a dry day, takes that month's
##' values pooled over the classes.
##'
##' \code{pulse_series_markov()} draws each year in turn: its total from the
##' quantiles of the annual totals, and so its class; the wet days from that
##' class's chain, month by month; and each wet day's amount from that month's
##' quantiles. It then scales the year's amounts to the drawn total. A model
##' year has 365 days and starts on 1 January, at \code{from}. Each day's
##' amount is a pulse at the middle of its day, or, with \code{days_per_pulse}
##' above 1, the amounts of that many days are one pulse at the middle of the
##' days. A pulse of zero is left out. The draws come from R's random number
##' stream, so \code{set.seed()} fixes the series.
##'
##' @title Series of resource pulses fitted to a daily record
##' @param date Dates of the record's days, as \code{Date} or something
##'   \code{as.Date()} takes.
##' @param amount Amount on each day, \code{NA} where the day is missing, in
##'   the units the pulse constructor takes (TF24: metres of water).
##' @param wet_threshold A day is wet when its amount is above this.
##' @param min_wet The fewest days a month and class needs for a fit of its
##'   own.
##' @return \code{pulse_series_fit()}: a list with the quantile levels
##'   \code{probs}, the quantiles of the annual totals \code{annual}, their
##'   \code{terciles}, the 12 x 3 (month by class) matrices \code{p_ww} and
##'   \code{p_wd} of the chance of a wet day after a wet and a dry day, the
##'   12 x 3 x \code{length(probs)} array \code{depth} of quantiles of the
##'   amounts on wet days, and the complete \code{years} used.
##'   \code{pulse_series_markov()}: a data frame with columns \code{time} and
##'   \code{amount}.
##' @seealso \code{\link{pulse_series_poisson}}, \code{\link{rainfall_pulse}}
##' @export
##' @examples
##' set.seed(1)
##' date <- seq(as.Date("1991-01-01"), as.Date("2020-12-31"), by = "day")
##' wet <- runif(length(date)) < ifelse(format(date, "%m") %in% c("12", "01", "02"), 0.5, 0.1)
##' amount <- ifelse(wet, rexp(length(date), rate = 1 / 0.01), 0)
##' fit <- pulse_series_fit(date, amount)
##' s <- pulse_series_markov(0, 100, fit, days_per_pulse = 7)
##' sum(s$amount) / 100          # about the record's mean annual total
pulse_series_fit <- function(date, amount, wet_threshold = 0, min_wet = 30) {
  date <- as.Date(date)
  if (length(date) == 0L || length(date) != length(amount)) {
    stop("`date` and `amount` must be non-empty and the same length")
  }
  if (anyNA(date) || anyDuplicated(date) > 0L) {
    stop("`date` must have no missing or repeated days")
  }
  if (any(!is.finite(amount) & !is.na(amount)) || any(amount < 0, na.rm = TRUE)) {
    stop("`amount` must be finite and non-negative, or NA where a day is missing")
  }
  o <- order(date)
  date <- date[o]
  amount <- amount[o]
  year <- as.integer(format(date, "%Y"))

  ## A complete year has an observation on every day of its calendar year,
  ## counted against 365 or 366, so a record that starts or stops mid-year does
  ## not look complete.
  obs <- tapply(!is.na(amount), year, sum)
  yrs <- as.integer(names(obs))
  need <- 365 + (yrs %% 4 == 0 & (yrs %% 100 != 0 | yrs %% 400 == 0))
  years <- yrs[obs == need]
  if (length(years) == 0L) {
    stop("the record has no complete year")
  }
  keep <- year %in% years
  amount <- amount[keep]
  year <- year[keep]
  month <- as.integer(format(date[keep], "%m"))
  wet <- amount > wet_threshold

  annual <- as.numeric(tapply(amount, year, sum))
  terciles <- stats::quantile(annual, c(1 / 3, 2 / 3), type = 8, names = FALSE)
  class <- (1L + (annual > terciles[1]) + (annual > terciles[2]))[match(year, years)]

  ## The days of a complete year are consecutive, so the day before is the row
  ## before; the first day of a year has none.
  n <- length(wet)
  prev <- c(NA, wet[-n])
  prev[c(TRUE, year[-1] != year[-n])] <- NA

  probs <- c(seq(0, 0.99, by = 0.01), 0.995, 0.999, 1)
  fit_cell <- function(cell) {
    after <- cell & !is.na(prev)
    d <- amount[cell & wet]
    list(n = c(length(d), sum(after & prev), sum(after & !prev)),
         p_ww = mean(wet[after & prev]),
         p_wd = mean(wet[after & !prev]),
         depth = if (length(d) > 0L) {
           stats::quantile(d, probs, type = 8, names = FALSE)
         } else {
           rep(0, length(probs))
         })
  }
  p_ww <- p_wd <- matrix(0, 12, 3)
  depth <- array(0, c(12, 3, length(probs)))
  for (m in seq_len(12)) {
    pooled <- fit_cell(month == m)
    for (k in seq_len(3)) {
      cell <- fit_cell(month == m & class == k)
      if (any(cell$n < min_wet)) {
        cell <- pooled
      }
      p_ww[m, k] <- cell$p_ww
      p_wd[m, k] <- cell$p_wd
      depth[m, k, ] <- cell$depth
    }
  }
  ## A month with no rain in the record leaves NaN: it is never wet.
  p_ww[!is.finite(p_ww)] <- 0
  p_wd[!is.finite(p_wd)] <- 0

  list(probs = probs,
       annual = stats::quantile(annual, probs, type = 8, names = FALSE),
       terciles = terciles, p_ww = p_ww, p_wd = p_wd, depth = depth,
       years = years)
}

##' @param fit A fit from \code{pulse_series_fit()}.
##' @param days_per_pulse Days whose amounts make one pulse: 1 for a pulse each
##'   wet day, 7 for one each week.
##' @rdname pulse_series_fit
##' @export
pulse_series_markov <- function(from, to, fit, days_per_pulse = 1) {
  check_pulse_span(from, to)
  if (!(length(days_per_pulse) == 1L && is.finite(days_per_pulse) &&
        days_per_pulse >= 1 && days_per_pulse == round(days_per_pulse))) {
    stop("`days_per_pulse` must be a whole number of days, at least 1")
  }
  ## Whole days only, so every pulse lies in (from, to].
  n_days <- floor((to - from) * 365 + 1e-9)
  if (n_days < 1) {
    stop("`from` and `to` must span at least one day")
  }
  depth <- unlist(lapply(seq_len(ceiling(n_days / 365)),
                         function(i) pulse_series_markov_year(fit)))
  depth <- depth[seq_len(n_days)]

  first <- seq(1, n_days, by = days_per_pulse)
  last <- pmin(first + days_per_pulse - 1, n_days)
  amount <- as.numeric(tapply(depth, (seq_len(n_days) - 1) %/% days_per_pulse, sum))
  time <- from + (first - 1 + last) / (2 * 365)
  keep <- amount > 0
  data.frame(time = time[keep], amount = amount[keep])
}

## One model year of daily amounts drawn from a fit, zero on dry days. It draws
## the year's total, then each day's chance of rain, then each day's amount.
pulse_series_markov_year <- function(fit) {
  u_annual <- stats::runif(1)
  u_occur <- stats::runif(365)
  u_depth <- stats::runif(365)
  quantile_at <- function(values, u) {
    stats::approx(fit$probs, values, xout = u, rule = 2)$y
  }
  target <- quantile_at(fit$annual, u_annual)
  k <- 1L + (target > fit$terciles[1]) + (target > fit$terciles[2])
  p_ww <- fit$p_ww[, k]
  p_wd <- fit$p_wd[, k]
  month <- rep(seq_len(12), c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))

  ## The first day has no day before it, so it is wet with the chain's
  ## stationary chance for January.
  wet <- logical(365)
  denom <- 1 + p_wd[1] - p_ww[1]
  wet[1] <- u_occur[1] < (if (denom > 0) p_wd[1] / denom else p_wd[1])
  for (d in 2:365) {
    wet[d] <- u_occur[d] < (if (wet[d - 1]) p_ww[month[d]] else p_wd[month[d]])
  }
  depth <- numeric(365)
  for (m in unique(month[wet])) {
    i <- which(wet & month == m)
    depth[i] <- quantile_at(fit$depth[m, k, ], u_depth[i])
  }

  ## The year is scaled to its own total and to nothing else: a run is a finite
  ## sample of years, and keeps the variation of one.
  if (sum(depth) > 0) {
    depth * (target / sum(depth))
  } else {
    depth[which.min(u_occur)] <- target
    depth
  }
}
