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
                      do.call(rbind, pad_matrix(x)))
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

  schedule <- schedule_or_null(schedule, p)

  control <- p$control$parameters
  eps <- control$schedule_eps
  progress <- control$schedule_progress

  history <- list()

  for (i in seq_len(control$schedule_nsteps)) {
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

##' Generate default cohort schedule.  This has times compressed
##' towards young patch ages, which usually need finer resolution than
##' later patche ages.  Currently there is no support for arbitrary
##' times: just build the schedule manually for now.
##'
##' @title Generate Default Cohort Schedule
##' @param p A Parameters object
##' @param t_max Maximum time to run schedule for.  If omitted (the
##' usual case) then time proceeds until there is a "very small"
##' chance that the patch is still extant, based on the control
##' parameter \code{schedule_default_patch_survival}.
##' @author Rich FitzJohn
##' @export
default_cohort_schedule <- function(p, t_max=NULL) {
  if (is.null(t_max)) {
    patch_survival <- p$control$parameters$schedule_default_patch_survival
    t_max <- p$disturbance$cdf(patch_survival)
  }
  times <- cohort_introduction_times(t_max)
  if (length(times) < 2) {
    stop("Did not generate at least two times, surprisingly")
  }
  n <- p$size
  sched <- new(CohortSchedule, n)
  sched$all_times <- rep(list(times[-(length(times))]), n)
  sched$max_time <- last(times)
  sched
}

## Internally used function.  This is a pattern that I've used before,
## but can't remember what I called it.  It might change in future so
## don't rely on it.  Basically: if we were given a schedule, copy
## it.  If we weren't, generate a default one from the parameters.
schedule_or_null <- function(schedule, p) {
  if (is.null(schedule)) {
    schedule <- default_cohort_schedule(p)
  } else {
    ## Don't modify what we're given:
    schedule <- schedule$copy()
  }
  schedule
}

##' Build a schedule of cohort times by iteratively finding the
##' position that causes the greatest change in fitness.  In contrast
##' with \code{\link{build_schedule}}, this is not really designed for
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
build_schedule_fitness <- function(p, schedule=NULL) {
  build_schedule_refine <- function(times, w, ebt, cores=1,
                                    verbose=FALSE) {
    run_with_insert <- function(i)
      run_with_times(insert_time(i, times), ebt)
    insert_time <- function(i, t) {
      j <- seq_len(i)
      c(t[j], (t[i] + t[i+1])/2, t[-j])
    }

    res <- unlist(mclapply(seq_len(length(times) - 1), run_with_insert,
                           mc.cores=cores))
    j <- which.max(abs(res - w))
    if (verbose)
      message(sprintf("Inserting time at %2.5f [%d]; fitness = %2.5f",
                      mean(times[j:(j+1)]), j, res[j]))

    list(times=insert_time(j, times), idx=j+1,
         seed_rain=cbind("in"=p$seed_rain, "out"=ebt$fitnesses))
  }
  run_with_times <- function(times, ebt) {
    ebt$reset()
    ebt$set_times(times, 1L)
    ebt$run()
    ebt$fitness(1L)
  }

  schedule <- schedule_or_null(schedule, p)

  control <- p$control$parameters
  eps <- control$schedule_eps
  progress <- control$schedule_progress

  ebt <- new(EBT, p)
  ebt$cohort_schedule <- schedule
  times <- ebt$cohort_schedule$times(1)

  w <- run_with_times(times, ebt)
  obj <- list(times=times,
              seed_rain=cbind("in"=p$seed_rain, "out"=ebt$fitnesses))
  history <- list(list(obj))

  for (i in seq_len(control$schedule_nsteps)) {
    obj <- build_schedule_refine(obj$times, obj$seed_rain[,"out"],
                                 ebt, verbose=verbose)
    history <- c(history, list(obj))
  }

  times <- obj$times
  attr(times, "seed_rain") <- obj$seed_rain
  if (progress)
    attr(times, "progress") <- history
  times
}
