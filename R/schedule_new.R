## Temporary file while I sort out new code....

##' Generate default cohort schedule.  This has times compressed
##' towards young patch ages, which usually need finer resolution than
##' later patche ages.  Currently there is no support for arbitrary
##' times: just build the schedule manually for now.
##'
##' @title Generate Default Cohort Schedule
##' @param p A Parameters object
##' @author Rich FitzJohn
##' @export
default_cohort_schedule <- function(p) {
  patch_survival <- p$control$parameters$schedule_default_patch_survival
  t_max <- p$disturbance$cdf(patch_survival)
  times <- cohort.introduction.times(t_max)
  if (length(times) < 2) {
    stop("Did not generate at least two times, surprisingly")
  }
  n <- p$size
  sched <- new(CohortSchedule, n)
  sched$all_times <- rep(list(times[-(length(times))]), n)
  sched$max_time <- last(times)
  sched
}

##' Build an appropriately refined schedule.
##'
##' There are control options (within the \code{Parameters} object)
##' that affect how this function runs, in particular
##' \code{schedule_nsteps} and \code{schedule_eps} control how refined
##' the schedule will end up, and \code{schedule_verbose} controls if
##' details are printed to the screen during construction.
##'
##' @title Build Cohort Schedule
##' @param p Parameters object
##' @param schedule An optional CohortSchedule, with a schedule to
##' start from.  If \code{NULL} (the default) a reasonable default
##' will be made using \code{\link{default_cohort_schedule}}
##' @return A CohortSchedule object
##' @author Rich FitzJohn
##' @export
build_schedule <- function(p, schedule=NULL) {
  run_with_lai <- function(schedule, p) {
    ebt <- new(EBT, p)
    ebt$cohort_schedule <- schedule

    err_lai <- rep(list(NULL), p$size)
    while (!ebt$complete) {
      added <- ebt$run_next()
      for (idx in added) {
        err_lai[[idx]] <-
          c(err_lai[[idx]], list(ebt$leaf_area_error(idx)))
      }
    }

    err_lai <- lapply(err_lai, function(x)
                      do.call(rbind, pad.matrix(x)))
    err_w <- lapply(seq_len(p$size), function(idx)
                    ebt$fitness_error(idx))
    f <- function(m)
      suppressWarnings(apply(m, 2, max, na.rm=TRUE))
    total <- lapply(seq_len(p$size), function(idx)
                    f(rbind(err_lai[[idx]], err_w[[idx]])))

    ## We'll often want the exact ode times that were used: this
    ## gathers them up, but does *not* set it up that they will be
    ## necessarily used.
    schedule$ode_times <- ebt$ode_times

    list(fitness=ebt$fitnesses,
         err=list(lai=err_lai, w=err_w, total=total))
  }
  split_times <- function(times, i) {
    ## Upwind splitting scheme only, which means that we will never
    ## split the last interval [assuming OK for now].  Inefficiently
    ## interleaved with sort().  These issues can change easily enough
    ## later.  The aim is making sure that we don't introduce the same
    ## point twice; one from upstream and one from downstream.
    dt <- diff(times)
    i <- which(i)
    sort(c(times, times[i] - dt[i-1]/2))
  }

  if (is.null(schedule)) {
    schedule <- default_cohort_schedule(p)
  } else {
    ## Don't modify what we're given:
    schedule <- schedule$copy()
  }

  control <- p$control$parameters
  nsteps <- control$schedule_nsteps
  eps <- control$schedule_eps
  progress <- control$schedule_progress

  history <- list()

  for (i in seq_len(nsteps)) {
    res <- run_with_lai(schedule, p)
    seed_rain <- cbind("in"=p$seed_rain, "out"=res[["fitness"]])
    split <- lapply(res$err$total, function(x) x > eps)

    ## TODO: I'm really not sure why we'd save the history.  There is
    ## probably a reason where this is used, but not sure what it is!
    ## If I don't find it soon, expect this code/functionality to be
    ## removed.
    if (progress) {
      history <- c(history, list(list(schedule=schedule$copy(),
                                      split=lapply(split, which),
                                      err=res$err,
                                      seed_rain=seed_rain)))
    }

    ## Prepare for the next iteration:
    for (idx in seq_len(p$size)) {
      schedule$set_times(split_times(schedule$times(idx), split[[idx]]),
                         idx)
    }

    if (!any(unlist(split), na.rm=TRUE)) {
      break
    }
    if (control$schedule_verbose) {
      message(sprintf("%d: Splitting {%s} times (%s)",
                      i,
                      paste(sapply(split, sum),    collapse=","),
                      paste(sapply(split, length), collapse=",")))
    }
  }

  if (progress) {
    attr(schedule, "progress") <- history
  }
  attr(schedule, "seed_rain") <- seed_rain

  schedule
}
