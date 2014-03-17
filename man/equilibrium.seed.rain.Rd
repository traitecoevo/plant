\name{equilibrium.seed.rain}
\alias{equilibrium.seed.rain}
\title{Determine Equilibrium Seed Rain}
\usage{
equilibrium.seed.rain(p, schedule, nsteps, build.args = list(),
  large.seed_rain.change = 10, progress = FALSE, verbose = TRUE)
}
\arguments{
  \item{p}{Parameters object}

  \item{schedule}{Schedule to run EBT with.  Will be
  modified with the final set of times.}

  \item{nsteps}{Number of iterations of the algorithm}

  \item{build.args}{Arguments to be passed through to
  \code{\link{build.schedule}}, as a list.}

  \item{large.seed_rain.change}{A change in seed rain
  "large enough" to trigger recalculation of the cohort
  schedule.  If the change is less than this value then we
  start the schedule search from the previous schedule.}

  \item{progress}{Return information on the convergence
  along with the output?}

  \item{verbose}{Be verbose about how the search is going?}
}
\value{
A list with elements \code{seed.in} (input seed rain on
last iteration), \code{seed.out} (output seed rain on last
iteration) and \code{times} (vector of cohort times, with
the final element being the stopping time).
}
\description{
Determine the equilibrium seed rain for a single strategy
}
\details{
Currently set up for just a single strategy; mostly because
build.schedule() only works with a single strategy, but
less importantly because the printing will not work
gracefully with multiple strategy.  None of these are
difficult to deal with though.

This does not use any clever root finding algorithm to try
and speed up convergence, but simply iterates the solver.
Something like the secant method would speed this up
greatly and help in divergent cases.
}
\author{
Rich FitzJohn
}

