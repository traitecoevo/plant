##' Determine the equilibrium seed rain for a single strategy
##'
##' Currently set up for just a single strategy; mostly because
##' build.schedule() only works with a single strategy, but less
##' importantly because the printing will not work gracefully with
##' multiple strategy.  None of these are difficult to deal with
##' though.
##'
##' This does not use any clever root finding algorithm to try and
##' speed up convergence, but simply iterates the solver.  Something
##' like the secant method would speed this up greatly and help in
##' divergent cases.
##'
##' @title Determine Equilibrium Seed Rain
##' @param p Parameters object
##' @param nsteps Number of iterations of the algorithm
##' @param times Vector of times to seed the cohort schedule with
##' (they will be computed with \code{\link{build.schedule}}
##' @param build.args Arguments to be passed through to
##' \code{\link{build.schedule}}, as a list.
##' @param large.seed_rain.change A change in seed rain "large enough" to
##' trigger recalculation of the cohort schedule.  If the change is
##' less than this value then we start the schedule search from the
##' previous schedule.
##' @param progress Return information on the convergence along with
##' the output?
##' @param verbose Be verbose about how the search is going?
##' @return A list with elements \code{seed.in} (input seed rain on
##' last iteration), \code{seed.out} (output seed rain on last
##' iteration) and \code{times} (vector of cohort times, with the
##' final element being the stopping time).
##' @author Rich FitzJohn
##' @export
equilibrium.seed.rain <- function(p, times, nsteps, build.args=list(),
                                  large.seed_rain.change=10,
                                  progress=FALSE, verbose=TRUE) {
  p <- p$clone() # don't modify what we're given
  times.default <- times
  build.args <- modifyList(list(nsteps=20, eps=1e-3, verbose=TRUE),
                           build.args)
  build <- function(times)
    build.schedule(p, times, build.args$nsteps, build.args$eps,
                   progress=FALSE, verbose=build.args$verbose)

  ## NOTE: The schedule of times here only gets larger.
  history <- list()
  for (i in seq_len(nsteps)) {
    times <- build(times)
    res <- list(seed_rain = attr(times, "seed_rain", exact=TRUE),
                times     = as.numeric(times))
    history <- c(history, list(res))
    seed_rain <- res[["seed_rain"]]
    change <- seed_rain[,"out"] - seed_rain[,"in"]

    p$seed_rain <- seed_rain[,"out"]
    if (abs(change) > large.seed_rain.change)
      times <- times.default

    if (verbose)
      message(sprintf("*** %d: %2.5f -> %2.5f (delta = %2.5f)", i,
                      seed_rain[,"in"], seed_rain[,"out"], change))
  }

  if (progress)
    attr(res, "progress") <- history
  res
}
