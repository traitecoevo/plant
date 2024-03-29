% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/k93.R
\name{make_K93_hyperpar}
\alias{make_K93_hyperpar}
\title{Hyperparameters for K93 physiological model}
\usage{
make_K93_hyperpar(
  b_0 = 0.059,
  b_1 = 0.012,
  b_2 = 0.00041,
  c_0 = 0.008,
  c_1 = 0.00044,
  d_0 = 0.00073,
  d_1 = 0.044,
  eta = 12,
  k_I = 0.01
)
}
\arguments{
\item{b_0}{Growth intercept year-1}

\item{b_1}{Growth asymptote year-1.(ln cm)-1}

\item{b_2}{Growth suppression rate m2.cm-2.year-1}

\item{c_0}{Mortality intercept year-1}

\item{c_1}{Mortality suppression rate m2.cm-2.year-1}

\item{d_0}{Recruitment rate (cm2.year-1)}

\item{d_1}{Recruitment suppression rate (m2.cm-2)}

\item{eta}{Crown shape parameter}

\item{k_I}{Extinction coefficient used when estimating competitive effect}
}
\description{
Construct hyperparameter object for K93 physiological model
}
