% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{run_plant_benchmarks}
\alias{run_plant_benchmarks}
\title{Check performance on current system using package bench}
\usage{
run_plant_benchmarks(
  strategy_types = list(FF16 = FF16_Strategy, FF16w = FF16w_Strategy, K93 = K93_Strategy),
  iterations = 1
)
}
\arguments{
\item{strategy_types}{A list of name strategy types to be tests}

\item{iterations}{The number of iterations to be run}
}
\value{
A dataframe of results
}
\description{
Check performance on current system using package bench
}
