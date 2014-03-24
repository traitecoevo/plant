\name{landscape}
\alias{landscape}
\alias{landscape.empty}
\title{Fitness Landscape}
\usage{
landscape(lma, p, schedule)

landscape.empty(lma, p, schedule)
}
\arguments{
  \item{lma}{Vector of lma values to compute mutant
  fitness}

  \item{p}{Parameters object.  Needs to contain residents
  with their incoming seed rain.}

  \item{schedule}{Schedule of times for the residents.}
}
\value{
Vector with the output seed rain.  Mutants have an
arbitrary seed rain of zero, so this is the rate of seed
production per capita.
}
\description{
Construct a fitness landscape.  Currently only works in the
LMA dimension, and with all other parameters set to their
defaults, and at predefined LMA values only rather than
adaptively (i.e., not actually useful yet).
}
\details{
Will change...
}
\author{
Rich FitzJohn
}

