\name{run.reference}
\alias{locate.evolve}
\alias{run.reference}
\title{Run Reference Model}
\usage{
  run.reference(path, p = NULL, path.evolve = NULL,
    verbose = TRUE)

  locate.evolve()
}
\arguments{
  \item{path}{Path to directory containing parameter files
  (see \code{\link{load.reference.output}} and
  \code{\link{load.reference.parameters}}}

  \item{p}{Optional \code{\link{Parameters}} object that
  will be written to \code{path} before running model
  (equivalent to "run model on these parameters in this
  directory")}

  \item{path.evolve}{Path to the "evolve" program.  If not
  specified, "evolve" is assumed to the in the
  \code{PATH}.}

  \item{verbose}{Logical indicating if the model output
  should be printed (TRUE by default).}
}
\description{
  Run the reference model, falster-traitdiversity
}
\details{
  This function is known to work with version d649534 of
  falster-traitdiversity (and will not work with any
  earlier version).
}
\author{
  Rich FitzJohn
}

