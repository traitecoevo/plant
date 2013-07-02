## EBT support functions.

##' Generate a suitable set of default cohort introduction times,
##' biased so that introductions are more closely packed at the
##' beginning of time, become increasingly spread out.
##'
##' See issue #30 for more details.
##'
##' @title Generate Default Cohort Introduction Times
##' @param max.time Time to generate introduction times up to (the
##' last introduction time will be at least \code{max.time}).
##' @param multiplier The rate of increase of step size with time.
##' The greater the number the faster step size will increase.
##' @param min.step.size The smallest gap between introduction times
##' (must be greater than zero, and will be the first introduction
##' time).
##' @param max.step.size The largest gap between introduction times
##' (may be infinite).
##' @return Vector of introduction times.
##' @export
##' @author Rich FitzJohn, adapted from original C++ code by Daniel
##' S. Falster.
cohort.introduction.times <- function(max.time, multiplier=0.2,
                                      min.step.size=1e-5,
                                      max.step.size=2.0) {
  if (min.step.size <= 0)
    stop("The minimum step size must be greater than zero")
  times <- numeric(0)
  dt <- time <- 0
  while (time <= max.time) {
    times <- c(times, time)
    dt <- 2^floor(log2(time * multiplier))
    time <- time + max(min(dt, max.step.size), min.step.size)
  }
  times
}
