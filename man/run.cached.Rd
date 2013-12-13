\name{run.cached}
\alias{run.cached}
\title{Run and Cache an Expression}
\usage{
  run.cached(expr, filename, regenerate = FALSE)
}
\arguments{
  \item{expr}{An expression to be evaluated}

  \item{filename}{A filename to store the result (will be
  saved via \code{\link{saveRDS}}}

  \item{regenerate}{Logical value indicating if results
  should be regenerated.}
}
\value{
  The result of evaluating \code{expr}, possibly loaded
  from \code{filename} rather than being rerun.
}
\description{
  Runs an expression and caches the result in a file.  Very
  basic, no checks.
}
\author{
  Rich FitzJohn
}

