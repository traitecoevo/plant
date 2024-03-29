% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logging.R
\name{plant_log_console}
\alias{plant_log_console}
\title{Activate logging with loggr}
\usage{
plant_log_console(
  file_name = "console",
  .message = FALSE,
  .warning = FALSE,
  .error = FALSE,
  ...
)
}
\arguments{
\item{file_name}{File to save output (default = console)}

\item{.message, .warning, .error}{Include messages, warnings or
errors?  By default (and in contrast to \code{loggr::log_file}
these are disabled here.}

\item{...}{Additional parameters passed to \code{loggr::log_file},
but \emph{not} \code{file_name} which is hard coded here to
\code{"console"}.}
}
\description{
Activate logging with loggr
}
\details{
By default plant prints little information about its progress.
This can be modified by enabling logging.  A formatter that is
different to the default \code{loggr::log_file} formatter is
selected here; it will print additional information that plant's
internal logging functions record.

"Schedule" events (splitting) are sent to the DEBUG stream,
everything else is sent to INFO.  All events have a "routine"
field added to them, which is useful if sent to a Redis server
(using \code{loggr.redis}).
}
