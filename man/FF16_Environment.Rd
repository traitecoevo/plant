% Generated by roxygen2: do not edit by hand
% Please edit documentation in r/RcppR6.R, r/ff16.R, r/ff16r.R
\name{FF16_Environment}
\alias{FF16_Environment}
\alias{FF16_make_environment}
\alias{FF16_fixed_environment}
\alias{FF16r_make_environment}
\title{FF16_Environment object}
\usage{
FF16_Environment(disturbance_mean_interval, seed_rain, control)

FF16_make_environment(p)

FF16_fixed_environment(e = 1, p = FF16_Parameters(), height_max = 150)

FF16r_make_environment(p)
}
\arguments{
\item{disturbance_mean_interval}{mean disturbance interval in years}

\item{seed_rain}{vector of seed rains}

\item{control}{A Control object}

\item{p}{A Parameters object}

\item{height_max}{= 150.0 maximum possible height in environment}

\item{e=1.0}{Value of environment}
}
\description{
FF16_Environment object

Construct a fixed environment for FF16 strategy
}
