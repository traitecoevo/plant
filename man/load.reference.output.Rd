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
  point. Each matrix has columns \code{h} (height) and
  \code{canopy.openness}.  These are always 50 rows long.}
  \item{popn}{A list of lists, one per strategy.  Each
  strategy's list has components: \describe{ \item{age}{Not
  sure!} \item{coh_m}{Not sure!} \item{coh_n}{Not sure!}
  \item{coh_r}{Not sure!} \item{bound_m}{Not sure!}
  \item{bound_s}{Not sure!} \item{bound_n}{Not sure!}
  \item{bound_r}{Not sure!} \item{d_bound_m}{Not sure!}
  \item{d_bound_s}{Not sure!} \item{d_bound_r}{Not sure!}
  \item{popn}{Not sure!} }} \item{parameters}{A
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

