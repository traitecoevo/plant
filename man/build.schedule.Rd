\name{build.schedule}
\alias{build.schedule}
\title{Build Cohort Schedule}
\usage{
  build.schedule(p, times, nsteps, eps, progress = FALSE,
    verbose = FALSE)
}
\arguments{
  \item{p}{Parameters object}

  \item{times}{Vector of times to start from}

  \item{nsteps}{Number of rounds of refinements to go
  through}

  \item{eps}{Numerical tolerance indicating where cohorts
  are refined enough.}

  \item{progress}{Save progress with the returned times?}

  \item{verbose}{Print information about progress as we
  go?}
}
\value{
  A vector of times
}
\description{
  Build an appropriately refined schedule for a single
  species.
}
\author{
  Rich FitzJohn
}

