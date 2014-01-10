## Code for building the cohort schedule by a series of refinements.
## This will change in future.  It also really only works for a single
## species at present.

##' Build an appropriately refined schedule for a single species.
##'
##' @title Build Cohort Schedule
##' @param p Parameters object
##' @param times Vector of times to start from
##' @param nsteps Number of rounds of refinements to go through
##' @param eps Numerical tolerance indicating where cohorts are
##' refined enough.
##' @param progress Save progress with the returned times?
##' @param verbose Print information about progress as we go?
##' @return A vector of times
##' @author Rich FitzJohn
##' @export
build.schedule <- function(p, times, nsteps, eps,
                           progress=FALSE, verbose=FALSE) {
  ## This gets the leaf area error out at each cohort introduction;
  ## that's the error criterion we care about.
  run.with.lai <- function(times, ebt) {
    idx <- 1
    ebt$reset()
    ebt$set_times(times, idx)
    err.lai <- list()
    while (ebt$cohort_schedule$remaining > 0) {
      ebt$run_next()
      err.lai <- c(err.lai, list(ebt$leaf_area_error(idx)))
    }
    err.lai <- do.call(rbind, pad.matrix(err.lai))
    err.w <- ebt$fitness_error(idx)
    total <- suppressWarnings(apply(rbind(err.lai, err.w), 2, max, na.rm=TRUE))

    list(lai=err.lai, w=err.w, total=total)
  }

  max.time <- last(times)
  ebt <- new(EBT, p)
  ebt$cohort_schedule <- schedule.from.times(times)
  times.intro <- ebt$cohort_schedule$times(1) # max.time dropped

  history <- list()
  for (i in seq_len(nsteps)) {
    err <- run.with.lai(times.intro, ebt)
    split <- err$total > eps

    if (progress)
      history <- c(history, list(list(times=c(times.intro, max.time),
                                      split=which(split), err=err)))

    ## Prepare for the next iteration:
    times.intro <- split.times(times.intro, split)
    times       <- c(times.intro, max.time)
    attr(times, "seed_rain") <- cbind("in" =p$seed_rain,
                                      "out"=ebt$fitnesses)

    if (!any(split, na.rm=TRUE))
      break
    if (verbose)
      message(sprintf("%d: Splitting %d times (%d)",
                      i, sum(split), length(times)))
  }

  if (progress)
    attr(times, "progress") <- history

  times
}

##' Build a schedule of cohort times by iteratively finding the
##' position that causes the greatest change in fitness.  In contrast
##' with \code{\link{build.schedule}}, this is not really designed for
##' use in real analyses, but instead to confirm that the schedules
##' work as expected.  It is \emph{extremely slow} because every step
##' requires introducing a cohort between every pair of cohorts.
##'
##' @title Build Schedule From Fitness Changes
##' @param p Parameters object
##' @param times Vector of times to start from
##' @param nsteps Number of new cohorts to add
##' @param progress Save progress with the returned times?
##' @param verbose Print information about progress as we go?
##' @return Numeric vector of times
##' @author Rich FitzJohn
##' @export
build.schedule.fitness <- function(p, times, nsteps,
                                   progress=FALSE, verbose=FALSE) {
  build.schedule.refine <- function(times, w, ebt, cores=1,
                                    verbose=FALSE) {
    run.with.insert <- function(i)
      run.with.times(insert.time(i, times), ebt)
    insert.time <- function(i, t) {
      j <- seq_len(i)
      c(t[j], (t[i] + t[i+1])/2, t[-j])
    }

    res <- unlist(mclapply(seq_len(length(times) - 1), run.with.insert,
                           mc.cores=cores))
    j <- which.max(abs(res - w))
    if (verbose)
      message(sprintf("Inserting time at %2.5f [%d]; fitness = %2.5f",
                      mean(times[j:(j+1)]), j, res[j]))

    list(times=insert.time(j, times), idx=j+1,
         seed_rain=cbind("in"=p$seed_rain, "out"=ebt$fitnesses))
  }
  run.with.times <- function(times, ebt) {
    ebt$reset()
    ebt$set_times(times, 1L)
    ebt$run()
    ebt$fitness(1L)
  }

  ebt <- new(EBT, p)
  ebt$cohort_schedule <- schedule.from.times(times)
  times <- ebt$cohort_schedule$times(1)

  w <- run.with.times(times, ebt)
  obj <- list(times=times,
              seed_rain=cbind("in"=p$seed_rain, "out"=ebt$fitnesses))
  history <- list(list(obj))

  for (i in seq_len(nsteps)) {
    obj <- build.schedule.refine(obj$times, obj$seed_rain[,"out"],
                                 ebt, verbose=verbose)
    history <- c(history, list(obj))
  }

  times <- obj$times
  attr(times, "seed_rain") <- obj$seed_rain
  if (progress)
    attr(times, "progress") <- history
  times
}

split.times <- function(times, i) {
  ## Upwind splitting scheme only, which means that we will never
  ## split the last interval [assuming OK for now].  Inefficiently
  ## interleaved with sort().  These issues can change easily enough
  ## later.  The aim is making sure that we don't introduce the same
  ## point twice; one from upstream and one from downstream.
  dt <- diff(times)
  i <- which(i)
  sort(c(times, times[i] - dt[i-1]/2))
}
