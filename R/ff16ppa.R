# Built from  R/ff16.R on Thu Oct 29 11:14:41 2020 using the scaffolder, from the strategy:  FF16
## We can probably actually do better than this with an S3 method on
## the actual strategy?  That would need to be organised by the
## templating though and that's stretched to the limit.

##' Create a FF16ppa Plant or Cohort
##' @title Create a FF16ppa Plant or Cohort
##' @param s A \code{\link{FF16ppa_Strategy}} object
##' @export
##' @rdname FF16ppa
##' @examples
##' pl <- FF16ppa_Individual()
##' pl$height
FF16ppa_Individual <- function(s=FF16ppa_Strategy()) {
  Individual("FF16ppa", "FF16_Env")(s)
}

#' Compute the whole plant light compensation point for a single
#' plant with FF16ppa strategy. Called via general function in plant.R
##' @export
##' @rdname FF16ppa
`lcp_whole_plant.Plant<FF16ppa>` <- function(p, ...) {
  FF16ppa_lcp_whole_plant(p, ...)
}

##' @export
##' @rdname FF16ppa
FF16ppa_Cohort <- function(s=FF16ppa_Strategy()) {
  Cohort("FF16ppa", "FF16_Env")(s)
}

##' @export
##' @rdname FF16ppa
FF16ppa_Species <- function(s=FF16ppa_Strategy()) {
  Species("FF16ppa", "FF16_Env")(s)
}

##' @export
##' @rdname FF16ppa
##' @param ... Arguments!
FF16ppa_Parameters <- function() {
  Parameters("FF16ppa","FF16_Env")()
}

##' @export
##' @rdname FF16ppa
##' @param p A \code{Parameters<FF16ppa,FF16_Env>} object
FF16ppa_Patch <- function(p) {
  Patch("FF16ppa", "FF16_Env")(p)
}

##' @export
##' @rdname FF16ppa
FF16ppa_SCM <- function(p) {
  SCM("FF16ppa", "FF16_Env")(p)
}

##' @export
##' @rdname FF16ppa
FF16ppa_StochasticSpecies <- function(s=FF16ppa_Strategy()) {
  StochasticSpecies("FF16ppa", "FF16_Env")(s)
}

##' @export
##' @rdname FF16ppa
FF16ppa_StochasticPatch <- function(p) {
  StochasticPatch("FF16ppa", "FF16_Env")(p)
}

##' @export
##' @rdname FF16ppa
FF16ppa_StochasticPatchRunner <- function(p) {
  StochasticPatchRunner("FF16ppa", "FF16_Env")(p)
}


## Helper:
##' @export
##' @rdname FF16_Environment
##' @param p A Parameters object
FF16ppa_make_environment <- function(p) {
  FF16_Environment(p$disturbance_mean_interval, p$seed_rain, p$k_I, p$control)
}

##' Construct a fixed environment for FF16ppa strategy
##'
##' @param e=1.0 Value of environment
##' @param p A Parameters object
##' @param height_max = 150.0 maximum possible height in environment
##' @rdname FF16_Environment
##'
##' @export
FF16ppa_fixed_environment <- function(e=1.0, p = FF16ppa_Parameters(), height_max = 150.0) {
  env <- FF16ppa_make_environment(p)
  env$set_fixed_environment(e, height_max)
  env
}


## This makes a pretend light environment over the plant height,
## slightly concave up, whatever.
FF16ppa_test_environment <- function(height, n=101, light_env=NULL,
                             n_strategies=1, seed_rain=0) {
  if (length(seed_rain) == 1) {
    seed_rain <- rep(seed_rain, length.out=n_strategies)
  }
  hh <- seq(0, height, length.out=n)
  if (is.null(light_env)) {
    light_env <- function(x) {
      exp(x/(height*2)) - 1 + (1 - (exp(.5) - 1))/2
    }
  }
  ee <- light_env(hh)
  interpolator <- Interpolator()
  interpolator$init(hh, ee)

  parameters <- FF16ppa_Parameters()
  parameters$strategies <- rep(list(FF16ppa_Strategy()), n_strategies)
  parameters$seed_rain <- seed_rain
  parameters$is_resident <- rep(TRUE, n_strategies)

  ret <- FF16ppa_make_environment(parameters)
  ret$canopy$canopy_interpolator <- interpolator
  attr(ret, "light_env") <- light_env
  ret
}


##' Hyperparameter function for FF16ppa physiological model
##' @title Hyperparameter function for FF16ppa physiological model
##' @param m A matrix of trait values, as returned by \code{trait_matrix}
##' @param s A strategy object
##' @param filter A flag indicating whether to filter columns. If TRUE, any numbers
##' that are within eps of the default strategy are not replaced.
##' @export
FF16ppa_hyperpar <- make_FF16_hyperpar()
