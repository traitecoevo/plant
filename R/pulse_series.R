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
  time <- from + seq_len(n) * interval
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
