\name{cohort.introduction.times}
\alias{cohort.introduction.times}
\title{Generate Default Cohort Introduction Times}
\usage{
  cohort.introduction.times(max.time, multiplier = 0.2,
    min.step.size = 1e-05, max.step.size = 2)
}
\arguments{
  \item{max.time}{Time to generate introduction times up to
  (the last introduction time will be at least
  \code{max.time}).}

  \item{multiplier}{The rate of increase of step size with
  time. The greater the number the faster step size will
  increase.}

  \item{min.step.size}{The smallest gap between
  introduction times (must be greater than zero, and will
  be the first introduction time).}

  \item{max.step.size}{The largest gap between introduction
  times (may be infinite).}
}
\value{
  Vector of introduction times.
}
\description{
  Generate a suitable set of default cohort introduction
  times, biased so that introductions are more closely
  packed at the beginning of time, become increasingly
  spread out.
}
\details{
  The reason for the stepped distribution is to keep step
  sizes as series of doublings.  Doing this limits the
  range of possible introduction times from an infinite set
  of possible values to a very limited subset of values
  (based on combinations of 1, 0.5, 0.25, 0.125 etc).  The
  reason for doing this is to minimise the number of unique
  introduction times across all species. The ODE stepper
  needs to stop at each point where a cohort is introduced.
  If each species was selecting a bunch of points that was
  essentially unique (compared to those selected for all
  other species), the number of unique cohort introductions
  times could get very large, requiring more ODE steps.
}
\author{
  Rich FitzJohn, adapted from original C++ code by Daniel
  S. Falster.
}

