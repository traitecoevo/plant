\name{build.schedule}
\alias{build.schedule}
\title{Build Cohort Schedule}
\usage{
  build.schedule(p, schedule, nsteps, eps,
    progress = FALSE, verbose = FALSE)
}
\arguments{
  \item{p}{Parameters object}

  \item{schedule}{A CohortSchedule, with a schedule to
  start from}

  \item{nsteps}{Number of rounds of refinements to go
  through}

  \item{eps}{Numerical tolerance indicating where cohorts
  are refined enough.}

  \item{progress}{Save progress with the returned times?}

  \item{verbose}{Print information about progress as we
  go?}
}
\value{
  A CohortSchedule object
}
\description{
  Build an appropriately refined schedule.
}
\author{
  Rich FitzJohn
}

