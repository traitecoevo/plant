\name{build.schedule.fitness}
\alias{build.schedule.fitness}
\title{Build Schedule From Fitness Changes}
\usage{
  build.schedule.fitness(p, nsteps, times,
    progress = FALSE, verbose = FALSE)
}
\arguments{
  \item{p}{Parameters object}

  \item{nsteps}{Number of new cohorts to add}

  \item{times}{Vector of times to start from}

  \item{progress}{Save progress with the returned times?}

  \item{verbose}{Print information about progress as we
  go?}
}
\value{
  Numeric vector of times
}
\description{
  Build a schedule of cohort times by iteratively finding
  the position that causes the greatest change in fitness.
  In contrast with \code{\link{build.schedule}}, this is
  not really designed for use in real analyses, but instead
  to confirm that the schedules work as expected.  It is
  \emph{extremely slow} because every step requires
  introducing a cohort between every pair of cohorts.
}
\author{
  Rich FitzJohn
}

