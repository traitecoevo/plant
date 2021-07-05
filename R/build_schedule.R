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
##' @return A Parameters object, with schedule components set.  The
##' output offspring produced is also available as an attribute
##' \code{birth_rate}.
##' @author Rich FitzJohn
##' @export
build_schedule <- function(p, state=NULL) {
  p <- validate(p)

  n_spp <- length(p$strategies)
  if (n_spp == 0L || !any(p$is_resident)) {
    stop("Can't build a schedule with no residents")
  }
  control <- p$control
  eps <- control$schedule_eps

  for (i in seq_len(control$schedule_nsteps)) {
    res <- run_scm_error(p, state)
    net_reproduction_ratios <- res[["net_reproduction_ratios"]]
    split <- lapply(res$err$total, function(x) x > eps)

    if (!any(unlist(split), na.rm=TRUE)) {
      break
    }

    ## Prepare for the next iteration:
    times <- res$schedule
    
    # split cohorts
    for (idx in seq_len(n_spp)) {
      
      # by introduction time
      times[[idx]] <- split_times(times[[idx]], split[[idx]])
      
      # or by initial size
      if(!is.null(state)) {
        state$species[[idx]] <- split_state(state$species[[idx]], 
                                            times[[idx]], split[[idx]],
                                            res$min_heights[[idx]])
      }
    }
    
    # set schedule for next patch
    p$cohort_schedule_times <- times
    msg <- sprintf("%d: Splitting {%s} times (%s)",
                   i,
                   paste(sapply(split, sum),    collapse=","),
                   paste(sapply(split, length), collapse=","))
    plant_log_debug(msg, routine="schedule", event="split", round=i)
  }

  p$cohort_schedule_ode_times <- res$ode_times
  ## Useful to record the last offspring produced:
  attr(p, "net_reproduction_ratios") <- net_reproduction_ratios

  return(list(parameters = p, state = state))
}

split_times <- function(times, i) {
  ## Upwind splitting scheme only, which means that we will never
  ## split the last interval [assuming OK for now].  Inefficiently
  ## interleaved with sort().  These issues can change easily enough
  ## later.  The aim is making sure that we don't introduce the same
  ## point twice; one from upstream and one from downstream.
  dt <- diff(times)
  i <- which(i)
  
  # can't split cohorts introduced in the same time step
  i <- i[dt[i] > 0]
  
  sort(c(times, times[i] - dt[i-1]/2))
}

split_state <- function(state, times, i, min_height) {
  ## oversimplified splitting scheme, introduce a new cohort 
  # half a step between two existing cohorts, even thought this leaves density
  # uncorrected and results in a different initial size distribution
  
  # can only split cohorts introduced at t0
  i <- which(i)
  
  # which is tricky when there's only one - we want to bring the 
  # first introduced (non-initialised) cohort up to meet it
  if(ncol(state) == 1 & i[1] == 2) 
    return(cbind(state, state - state / 2))
  
  # otherwise proceed normally
  i <- i[times[i] == 0]
  
  if(length(i) == 0)
    return(state)
  
  # append a blank cohort
  dh <- matrix(0, ncol = ncol(state), nrow = nrow(state), 
               dimnames = dimnames(state)) 
  
  dh["height", ] <- diff(c(state["height", ], min_height))
  
  st <- cbind(state, state[, i-1] + dh[, i-1] / 2)

  return(st[, order(st["height", ], decreasing = T)])
}

