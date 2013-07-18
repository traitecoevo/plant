\name{load.reference.output}
\alias{load.reference.output}
\title{Load Reference Model Output}
\usage{
  load.reference.output(path)
}
\arguments{
  \item{path}{The path containing the output; this is the
  directory that contains the files "Stoch.txt",
  "params.m", and the directory "res" that has the actual
  output files.}
}
\value{
  A list with components: \item{stand}{Not sure!}
  \item{patch_age}{\code{data.frame} with columns
  \code{age} (age in years), \code{density} (density of
  patches of that age), \code{gamma} (parameter of the
  disturbance function).  After that are columns
  \code{add.cohort.1}, \code{add.cohort.2}, etc, which are
  logical vectors indicating if a cohort for a strategy was
  added at that time point.} \item{seed_rain_out}{Vector of
  output seed rain, one entry per strategy/species}
  \item{light_env}{Light environment as one matrix per time
  point. Each matrix has columns \code{height} (height) and
  \code{canopy.openness}.  These are always 50 rows long.}
  \item{popn}{A list of lists, one per strategy.  Each
  strategy's list has components: \describe{ \item{age}{Age
  of the cohort at a given point in time.}
  \item{coh_m}{Leaf mass of the mean individual of a
  cohort} \item{coh_n}{Mean density (leaf mass, see eqn
  22)} \item{coh_r}{Mean seed production}
  \item{bound_m}{Leaf mass of the boundary individual of a
  cohort} \item{bound_s}{Boundary seed production}
  \item{bound_n}{Boundary density (leaf mass, see eqn 22)}
  \item{bound_r}{Boundary seed production}
  \item{d_bound_m}{Time derivative of \code{bound_m}}
  \item{d_bound_s}{Time derivative of \code{bound_s}}
  \item{d_bound_r}{time derivative of \code{bound_r}}
  \item{popn}{Matrix of physiological/demographic summary
  statistics, one row per time.} }} \item{parameters}{A
  \code{\link{Parameters}} object}
}
\description{
  Load the output of running the reference model
  (falster-traitdiversity).
}
\details{
  There is still work to do in sanitising the output ready
  for comparison, so the return value will change.

  This function is known to work with version d649534 of
  falster-traitdiversity (and will not work with any
  earlier version).
}
\author{
  Rich FitzJohn
}

