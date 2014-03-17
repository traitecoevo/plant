\name{default.schedule}
\alias{default.schedule}
\title{Simple Linear Cohort Schedule}
\usage{
default.schedule(nt, max.t)
}
\arguments{
  \item{nt}{Number of cohorts to introduce}

  \item{max.t}{Time that the simulation stops at}
}
\value{
A \code{CohortSchedule} object
}
\description{
Generate a default schedule of cohort introduction times.
At present this is just linear from 0..max.t.
}
\details{
Note that this is basically useless for practical purposes,
but will merge with \code{\link{cohort.introduction.times}}
at some point.
}
\author{
Rich FitzJohn
}

