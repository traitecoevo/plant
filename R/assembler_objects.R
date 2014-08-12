## The simplest thing is a species.  It has up to three elements:
##
## traits: a (named) list of species traits
## seed_rain: current incoming seed rain
## cohort_schedule_times: (optional) a set of times for cohort
##   introductions
##' @export
species <- function(traits, seed_rain=1, cohort_schedule_times=NULL) {
  ## TODO: check traits are named.
  ret <- list(traits=traits,
              seed_rain=seed_rain,
              cohort_schedule_times=cohort_schedule_times)
  class(ret) <- "assembly_species"
  ret
}

## The next most complicated thing is a community.  At it's simplest
## this is just a list of species.  However we provide a number of
## methods for working with them:
##   size(): number of species in the community at present
##   traits(): matrix of traits across all species (one row per species)
##   seed_rain(): vector of seed rain across all species
##   to_parameters(): Generate Parameters object for community
##   to_schedule(): Generate CohortSchedule object for community
## and for modifying them:
##   set_seed_rain(): Set seed rain given a new vector
##   set_schedule_times(): Set schedule times from a CohortSchedule object
##   drop(): Drop some species, given a logical vetor
##   add_traits(): Add new species, given a matrix of new traits
## and for doing more things:
##   run(): Run the simulation with the current community, updating
##     see rain and cohort schedule times
##   make_landscape(): Using results from the last run, generate a
##     function that can compute 1D landscapes.
##
## One downside of R6 is that we can't (currently) add a method to an
## existing class bit-by-bit.  Getting around this with local
## functions.
.R6_community <- local({
  initialize <- function(p0, trait_names, ..., bounds=NULL,
                         seed_rain_initial=1e-3, sys0=NULL) {
    parameters <<- p0$copy()
    seed_rain_initial <<- seed_rain_initial
    if (!is.null(sys0)) {
      if (length(...) > 0) {
        stop("Don't provide extra things with sys0")
      }
    } else {
      sys0 <- list(...)
    }
    ## TODO: Check that all elements are species, and that they have
    ## the same set of variable traits.
    if (!all(sapply(sys0, inherits, "species"))) {
      stop("All elements must be species")
    }
    sys <<- sys0
    trait_names <<- trait_names
    if (length(trait_names) > 1) {
      stop("Doesn't yet support multiple traits")
    }
    if (is.null(bounds)) {
      n <- length(trait_names)
      b <- cbind(lower=rep(-Inf, n), upper=rep(Inf, n))
      rownames(b) <- trait_names
      bounds <<- bounds
    } else {
      if (nrow(bounds) != length(trait_names)) {
        stop("Incorrect size bounds")
      }
      if (ncol(bounds) != 2) {
        stop("Bounds must have two columns")
      }
      bounds <<- bounds
    }
  }
  size <- function() {
    length(sys)
  }
  traits <- function(as_matrix=FALSE) {
    collect("traits", sys, empty=NULL, loop=lapply, each=unlist,
            after=if (as_matrix) rbind_list else identity)
  }
  seed_rain <- function() {
    collect("seed_rain", sys, empty=numeric(0))
  }
  to_parameters <- function() {
    p <- parameters$copy()
    p$clear() # start empty

    ## Ideally we'd use expand_parameters, but update it to support a
    ## matrix of traits/values
    ##   expand_parameters(trait, value, p)
    traits <- traits(FALSE)
    for (i in traits) {
      s <- p$strategy_default$copy()
      s$set_parameters(i)
      p$add_strategy(s)
    }
    p$seed_rain <- seed_rain()
    p
  }
  to_schedule <- function(value, p=to_parameters()) {
    schedule <- default_cohort_schedule(p)
    for (i in seq_len(size())) {
      if (!is.null(sys[[i]]$cohort_schedule_times)) {
        schedule$set_times(sys[[i]]$cohort_schedule_times, i)
      }
    }
    schedule
  }
  set_seed_rain <- function(value) {
    n <- size()
    if (n != length(value)) {
      stop(sprintf("Invalid length: expected %d, recieved %d",
                   n, length(value)))
    }
    for (i in seq_len(n)) {
      sys[[i]]$seed_rain <<- value[[i]]
    }
  }
  set_schedule_times <- function(value) {
    if (!inherits(value, "Rcpp_CohortSchedule")) {
      stop("Expected CohortSchedule object")
    }
    n <- size()
    if (value$n_species != n) {
      stop(sprintf("Invalid length: expected %d, recieved %d",
                   n, value$n_species))
    }
    for (i in seq_len(n)) {
      sys[[i]]$cohort_schedule_times <<- value$times(i)
    }
  }
  set_viable_bounds <- function() {
    message("Computing viable bounds")
    if (nrow(bounds) != 1) {
      stop("This is not going to work with multiple traits yet")
    }
    bounds <<- viable_fitness(trait_names, to_parameters(),
                              bounds=base::drop(bounds))
  }
  drop <- function(which) {
    if (is.logical(which)) {
      if (length(which) != size()) {
        stop(sprintf("Invalid length: expected %d, recieved %d",
                     size(), length(which)))
      }
      keep <- !which
    } else { # should handle integers some time
      stop("Invalid index")
    }
    sys <<- sys[keep]
  }
  add_traits <- function(x) {
    assertthat::assert_that(is.matrix(x))
    if (size() > 0) {
      existing <- traits(TRUE)
      assertthat::assert_that(ncol(x) == ncol(existing))
      assertthat::assert_that(identical(colnames(x),
                                        trait_names))
    }
    n <- nrow(x)
    if (n > 0) {
      tmp <- vector("list", n)
      for (i in seq_len(nrow(x))) {
        tmp[[i]] <- species(as.list(x[i,]), seed_rain_initial)
      }
      sys <<- append(sys, tmp)
    }
  }
  run <- function(compute_schedule=TRUE) {
    if (size() > 0) {
      p <- to_parameters()
      ## Check for non-NA seed rain here, or you get very
      ## cryptic error.
      if(any(is.na(p$seed_rain)))
        stop("NA in seed rain")
      if (compute_schedule) {
        res <- build_schedule(p, to_schedule(p))
      } else {
        ## Bit of fakery here to generate ode times
        message("Recomputing ode times")
        ebt <- run_ebt(p, to_schedule(p))
        res <- ebt$cohort_schedule$copy()
        attr(res, "seed_rain") <- cbind(out=ebt$fitnesses)
      }
      seed_rain_out <- unname(attr(res, "seed_rain", exact=TRUE)[,"out"])
      set_seed_rain(seed_rain_out)
      set_schedule_times(res)
      last_p <<- p
      last_schedule <<- res
    } else {
      last_p <<- NULL
      last_schedule <<- NULL
    }
  }

  run_to_equilibrium <- function() {
    if (size() > 0) {
      p <- to_parameters()
      ## Check for non-NA seed rain here, or you get very
      ## cryptic error.
      if(any(is.na(p$seed_rain)))
        stop("NA in seed rain")

      res <- equilibrium_seed_rain(p)
      set_seed_rain(rowMeans(res$seed_rain))
      set_schedule_times(res$schedule)
      last_p <<- p
      last_schedule <<- res$schedule
    } else {
      last_p <<- NULL
      last_schedule <<- NULL
    }
  }
  make_landscape <- function() {
    if (!is.null(last_p) && !is.null(last_schedule)) {
      p <- last_p$copy()
      schedule <- last_schedule$copy()
      f <- function(x) {
        landscape(trait_names, x, p, schedule)
      }
    } else if (size() == 0) {
      p <- to_parameters()
      f <- function(x) {
        landscape_empty(trait_names, x, p)
      }
    } else {
      f <- NULL
    }
    f
  }
  store_sys_attribute <- function(value, name) {
    attr(sys, name) <<- value
  }
  R6::R6Class("community",
              public=list(
                initialize=initialize,
                get_sys=function() sys,
                size=size,
                drop=drop,
                trait_names=NULL,
                traits=traits,
                bounds=NULL,
                add_traits=add_traits,
                seed_rain=seed_rain,
                ## Functions generating useful things
                to_parameters=to_parameters,
                to_schedule=to_schedule,
                set_viable_bounds=set_viable_bounds,
                make_landscape=make_landscape,
                ## Fuctions that modify things:
                set_seed_rain=set_seed_rain,
                set_schedule_times=set_schedule_times,
                run=run,
                run_to_equilibrium=run_to_equilibrium,
                ## This is a hack for now:
                store_sys_attribute=store_sys_attribute),
              private=list(
                sys=NULL,
                parameters=NULL,
                seed_rain_initial=NA_real_,
                last_p=NULL,
                last_schedule=NULL))
})
##' @export
community <- function(...) {
  .R6_community$new(...)
}

## Then the highest level for now: the assembler:
.R6_assembler <- local({
  initialize <- function(community0, births_sys, deaths_sys,
                         filename=NULL) {
    community <<- community0
    births_sys <<- births_sys
    deaths_sys <<- deaths_sys
    history <<- list()
    filename <<- filename
    append()
  }
  deaths <- function() {
    deaths_sys(community)
  }
  births <- function() {
    to_add <- births_sys(community)
    if (nrow(to_add) > 0) {
      community$add_traits(to_add)
    }
  }
  step <- function(type="single") {
    f <- switch(type,
      single=step_single,
      to_equilibrium=step_to_equilibrium,
      stop("unknown type ", type))
    f()
  }
  step_single <- function() {
    message(sprintf("*** Assembler: step %d, (%d strategies)",
                    length(history), community$size()))
    deaths()
    births()
    community$run()
    append()
  }
  step_to_equilibrium <- function() {
    message(sprintf("*** Assembler: step %d, (%d strategies)",
                    length(history), community$size()))
    deaths()
    births()
    community$run_to_equilibrium()
    append()
  }
  append <- function() {
    history <<- c(history, list(community$get_sys()))
    if (!is.null(filename)) {
      ok <- try(saveRDS(history, filename))
      if (inherits(ok, "try-error")) {
        warning("History saving has failed",
                immediate.=TRUE, call.=FALSE)
      }
    }
  }
  run_nsteps <- function(n, type="single") {
    for(i in seq_len(n)) {
       step(type)
    }
  }

  R6::R6Class("assembler",
              public=list(
                initialize=initialize,
                deaths=deaths,
                births=births,
                step=step,
                step_single=step_single,
                step_to_equilibrium=step_to_equilibrium,
                run_nsteps=run_nsteps,
                get_community=function() community,
                get_history=function() history
                ),
              private=list(
                community=NULL,
                births_sys=NULL,
                deaths_sys=NULL,
                history=NULL,
                filename=NULL,
                append=append
                ))
})
##' @export
assembler <- function(...) {
  .R6_assembler$new(...)
}

## This is a hack for now.  To make this work more smoothly, we want
## to save some additional parameters (e.g., trait_names) out with the
## community.
restore_community <- function(x, p0, bounds=NULL, recompute=TRUE) {
  if (!all(sapply(x, inherits, "assembly_species"))) {
    stop("Not all elements are assembly_species")
  }
  if (length(x) == 0) {
    stop("Can't restore empty community")
  }
  trait_names <- names(x[[1]]$traits)
  if (!all(sapply(x[-1], function(el)
                  identical(names(el$traits), trait_names)))) {
    stop("Traits vary across species")
  }
  p <- p0$copy()
  p$clear()
  sys0 <- community(p, trait_names)
  sys0$private$sys <- x
  if (recompute) {
    sys0$run(FALSE)
  }
  ## TODO: Need to update bounds here, but that really wants to have
  ## the bounds checked.  Far better to write them out in the file at
  ## the same time, really!
  sys0
}
