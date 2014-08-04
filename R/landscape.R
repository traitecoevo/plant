##' Construct a fitness landscape.  Currently only works in one
##' dimension, and with all other parameters set to their defaults,
##' and at predefined values only rather than adaptively (i.e., not
##' actually useful yet).
##'
##' Will change...
##' @title Fitness Landscape
##' @param trait Name of the single trait to change
##' @param values Vector of trait values to compute mutant fitness
##' @param p Parameters object.  Needs to contain residents with their
##' incoming seed rain.
##' @param schedule Schedule of times for the residents.
##' @return Vector with the output seed rain.  Mutants have an
##' arbitrary seed rain of zero, so this is the rate of seed
##' production per capita.
##' @author Rich FitzJohn
##' @export
landscape <- function(trait, values, p, schedule=NULL) {
  if (is.null(schedule)) {
    schedule <- build_schedule(p)
  }

  ## We want to hit exactly the same times as the previous time this
  ## was run (with just the residents, establishing the schedule).
  ## This should usually be just fine
  schedule$use_ode_times <- TRUE
  p_with_mutants <- expand_parameters(trait, values, p)

  schedule_with_mutants <- expand_schedule(schedule, length(values))
  ebt_with_mutants <- run_ebt(p_with_mutants, schedule_with_mutants)
  w_with_mutants <- ebt_with_mutants$fitnesses
  w_with_mutants[-seq_len(p$size)]
}

## TODO: Need to support variant strategy here and above.
##' @rdname landscape
##' @export
landscape_empty <- function(trait, values, p, schedule=NULL) {
  p_empty <- p$copy()
  p_empty$clear()
  p_empty <- expand_parameters(trait, values, p_empty)

  if (is.null(schedule)) {
    schedule_empty <- default_cohort_schedule(p_empty)
  } else {
    ## Build an appropriate schedule:
    schedule_empty <- new(CohortSchedule, length(values))
    schedule_empty$max_time  <- schedule$max_time
    schedule_empty$all_times <-
      rep(list(unique(sort(unlist(schedule$all_times)))), p_empty$size)
    schedule_empty$ode_times <- schedule$ode_times
    ## We'll only use these if use_ode_times was set to TRUE in the
    ## original schedule...
    schedule_empty$use_ode_times <- schedule$use_ode_times
  }

  ebt_empty <- run_ebt(p_empty, schedule_empty)
  ebt_empty$fitnesses
}

##' Expand schedule to include mutants.  All mutants get the same
##' schedule, equal to all the unique times that any resident was
##' introduced.  This results in more work than is really needed, but
##' should be reasonable most of the time.
##'
##' May change...
##' @title Expand Cohort Schedule to Include Mutants
##' @param schedule A \code{CohortSchedule} object, set up for the
##' residents.
##' @param n_mutant Number of mutant schedules to generate
##' @author Rich FitzJohn
##' @export
expand_schedule <- function(schedule, n_mutant) {
  n_resident <- schedule$n_species
  ret <- new(CohortSchedule, n_resident + n_mutant)
  ret$max_time <- schedule$max_time

  ## Copy residents over:
  for (i in seq_len(n_resident))
    ret$set_times(schedule$times(i), i)

  ## Introduce mutants at all unique times:
  times_mutant <- unique(sort(unlist(schedule$all_times)))
  for (i in seq_len(n_mutant)) {
    ret$set_times(times_mutant, n_resident + i)
  }

  if (length(schedule$ode_times) == 0) {
    ## Basically, we want to step on the same times that the real EBT
    ## did.  To do that we need to get times out of the EBT.  If you
    ## have the ebt instance handy, just do this:
    ##   schedule$ode_times <- ebt$times
    ## If you ran build_schedule this will be set already.
    stop("Please add ode_times from your last ebt run")
  }
  ret$ode_times <- schedule$ode_times

  ret

}

##' Expand parameters to include mutants.  All mutants are built off
##' the default strategy.
##'
##' @title Expand Parameters to Include Mutants
##' @param trait Name of the single trait to change
##' @param values Vector of trait values for the mutants
##' @param p Parameters object
##' @param mutant Are these mutant strategies? (Default: yes they
##' are).
##' @author Rich FitzJohn
##' @export
expand_parameters <- function(trait, values, p, mutant=TRUE) {
  p <- p$copy()
  strategy <- p$strategy_default$copy()

  ## TODO: Generalise this out (easy)
  for (i in values) {
    new_strategy <- strategy$copy()
    new_strategy$set_parameters(structure(list(i), names=trait))
    if (mutant) {
      p$add_strategy_mutant(new_strategy)
    } else {
      p$add_strategy(new_strategy)
    }
  }

  p
}
