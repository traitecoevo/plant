\name{schedule.from.times}
\alias{schedule.from.times}
\title{Create CohortSchedule From Vector of Times}
\usage{
  schedule.from.times(times)
}
\arguments{
  \item{times}{A numeric vector of increasing times}
}
\value{
  A \code{\link{CohortSchedule}} object.
}
\description{
  Generate a CohortSchedule object from a series of times.
  The last time is ommited but will be set as the
  \code{$max_time}.
}
\details{
  Only works for a single species at the moment.
}
\author{
  Rich FitzJohn
}

