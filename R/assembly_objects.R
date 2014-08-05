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
  initialize <- function(p0, ..., seed_rain_initial=1e-3, sys0=NULL) {
    parameters <<- p0$copy()
    seed_rain_initial <<- seed_rain_initial
    if (!is.null(sys)) {
      if (length(...) > 0) {
        stop("Don't provide extra things with sys0")
      }
    } else {
      sys0 <- list(...)
    }
    ## TODO: Check that all elements are species, and that they have
    ## the same set of variable traits.
    if (!all(sapply(sys0, function(x) inherits, "species"))) {
      stop("All elements must be species")
    }
    sys <<- sys0
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
                                        colnames(existing)))
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
  run <- function() {
    if (size() > 0) {
      p <- to_parameters()
      ## TODO: Check for non-NA seed rain here, or you get very
      ## cryptic error.
      res <- build_schedule(p, to_schedule(p))
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
  make_landscape <- function() {
    if (is.null(last_p) || is.null(last_schedule)) {
      stop("Didn't find p/schedule from previous run")
    }
    p <- last_p$copy()
    schedule <- last_schedule$copy()
    trait <- colnames(traits(TRUE))
    if (length(trait) != 1) {
      stop("Landscape only works with one trait")
    }
    function(x) {
      landscape(trait, x, p, schedule)
    }
  }
  R6::R6Class("community",
              public=list(
                initialize=initialize,
                get_sys=function() sys,
                size=size,
                drop=drop,
                add_traits=add_traits,
                traits=traits,
                seed_rain=seed_rain,
                ## Functions generating useful things
                to_parameters=to_parameters,
                to_schedule=to_schedule,
                make_landscape=make_landscape,
                ## Fuctions that modify things:
                set_seed_rain=set_seed_rain,
                set_schedule_times=set_schedule_times,
                run=run),
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
##
##
.R6_assembler <- local({
  initialize <- function(community0, births_sys, deaths_sys) {
    community <<- community0
    births_sys <<- births_sys
    deaths_sys <<- deaths_sys
  }
  deaths <- function() {
    deaths_sys(community)
  }
  births <- function(must_grow=FALSE) {
    births_sys(community, must_grow)
  }
  step <- function() {
    deaths()
    births()
    ## This is where a pre-flight check of mutant fitness would go,
    ## easy if we have the fitness landscapes coming out of the
    ## previous step.
    community$run()
  }
  R6::R6Class("assembler",
              public=list(
                initialize=initialize,
                deaths=deaths,
                births=births,
                step=step,
                get_community=function() community
                ),
              private=list(
                community=NULL,
                births_sys=NULL,
                deaths_sys=NULL
                ))
})
##' @export
assembler <- function(...) {
  .R6_assembler$new(...)
}

## Then a specialised version of this based around some simple ideas
## of stochastic assembly.  This takes relatively few arguments --
## bounds is the only required one really.
.R6_assembler_stochastic <- local({
  initialize <- function(community0,
                         bounds, n_mutants=1L, n_immigrants=1L, vcv=NULL,
                         vcv_p=0.001,
                         seed_rain_eps=1e-3) {
    community <<- community0
    if (is.null(vcv)) {
      vcv <- vcv_p * diag(nrow(bounds)) * as.numeric(diff(t(log(bounds))))
    }
    births_sys <<- make_births(n_mutants, vcv, n_immigrants, bounds)
    deaths_sys <<- make_deaths(seed_rain_eps)
  }
  R6::R6Class("assembler_stochastic",
              inherit=.R6_assembler,
              public=list(initialize=initialize))
})
##' @export
assembler_stochastic <- function(...) {
  .R6_assembler_stochastic$new(...)
}
