\name{schedule.from.times}
\alias{schedule.from.times}
\title{Create CohortSchedule From Vector of Times}
\usage{
schedule.from.times(times, n = 1L)
}
\arguments{
  \item{times}{A numeric vector of increasing times}

  \item{n}{The number of species.}
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

