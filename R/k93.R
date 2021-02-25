# Built from  R/ff16.R on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16

## We can probably actually do better than this with an S3 method on
## the actual strategy?  That would need to be organised by the
## templating though and that's stretched to the limit.

##' Create a K93 Individual or Cohort
##' @title Create a K93 Individual or Cohort
##' @param s A \code{\link{K93_Strategy}} object
##' @export
##' @rdname K93
##' @examples
##' pl <- K93_Individual()
##' pl$height
K93_Individual <- function(s=K93_Strategy()) {
  Individual("K93", "K93_Env")(s)
}

#' Compute the whole plant light compensation point for a single
#' plant with K93 strategy. Called via general function in plant.R
##' @export
##' @rdname K93
`lcp_whole_plant.Individual<K93>` <- function(p, ...) {
  K93_lcp_whole_plant(p, ...)
}

##' @export
##' @rdname K93
K93_Cohort <- function(s=K93_Strategy()) {
  Cohort("K93", "K93_Env")(s)
}

##' @export
##' @rdname K93
K93_Species <- function(s=K93_Strategy()) {
  Species("K93", "K93_Env")(s)
}

##' @export
##' @rdname K93
##' @param ... Arguments!
K93_Parameters <- function() {
  Parameters("K93","K93_Env")()
}

##' @export
##' @rdname K93
##' @param p A \code{Parameters<K93,K93_Env>} object
K93_Patch <- function(p) {
  Patch("K93", "K93_Env")(p)
}

##' @export
##' @rdname K93
K93_SCM <- function(p) {
  SCM("K93", "K93_Env")(p)
}

##' @export
##' @rdname K93
K93_StochasticSpecies <- function(s=K93_Strategy()) {
  StochasticSpecies("K93", "K93_Env")(s)
}

##' @export
##' @rdname K93
K93_StochasticPatch <- function(p) {
  StochasticPatch("K93", "K93_Env")(p)
}

##' @export
##' @rdname K93
K93_StochasticPatchRunner <- function(p) {
  StochasticPatchRunner("K93", "K93_Env")(p)
}


## Helper:
##' @export
##' @rdname K93_Environment
##' @param p A Parameters object
K93_make_environment <- function(p) {
  K93_Environment(p$disturbance_mean_interval, p$seed_rain, p$k_I, p$control)
}

##' Construct a fixed environment for K93 strategy
##'
##' @param e=1.0 Value of environment 
##' @param p A Parameters object
##' @param height_max = 150.0 maximum possible height in environment
##' @rdname K93_Environment
##'
##' @export
K93_fixed_environment <- function(e=1.0, p = K93_Parameters(), height_max = 300.0) {
  env <- K93_make_environment(p)
  env$set_fixed_environment(e, height_max)
  env
}


## This makes a pretend light environment over the plant height,
## slightly concave up, whatever.
K93_test_environment <- function(height, n=101, light_env=NULL,
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

  parameters <- K93_Parameters()
  parameters$strategies <- rep(list(K93_Strategy()), n_strategies)
  parameters$seed_rain <- seed_rain
  parameters$is_resident <- rep(TRUE, n_strategies)

  ret <- K93_make_environment(parameters)
  ret$canopy$canopy_interpolator <- interpolator
  attr(ret, "light_env") <- light_env
  ret
}


##' Hyperparameters for K93 physiological model
##' @title Hyperparameters for K93 physiological model
##' @export
##' @rdname K93_hyperpar
make_K93_hyperpar <- function(
        b_0 = 0.059,    # Growth intercept year-1
        b_1 = 0.012,    # Growth asymptote year-1.(ln cm)-1
        b_2 = 0.00041,  # Growth suppression rate m2.cm-2.year-1
        c_0 = 0.008,    # Mortality intercept year-1
        c_1 = 0.00044,  # Mortality suppression rate m2.cm-2.year-1
        d_0 = 0.00073,  # Recruitment rate (cm2.year-1)
        d_1 = 0.044    # Recruitment suppression rate (m2.cm-2)
  ) {
  assert_scalar <- function(x, name=deparse(substitute(x))) {
    if (length(x) != 1L) {
      stop(sprintf("%s must be a scalar", name), call. = FALSE)
    }
  }

  assert_scalar(b_0)
  assert_scalar(b_1)
  assert_scalar(b_2)
  assert_scalar(c_0)
  assert_scalar(c_1)
  assert_scalar(d_0)
  assert_scalar(d_1)

  function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }
    
    m
  }

}


##' Hyperparameter function for K93 physiological model
##' @title Hyperparameter function for K93 physiological model
##' @param m A matrix of trait values, as returned by \code{trait_matrix}
##' @param s A strategy object
##' @param filter A flag indicating whether to filter columns. If TRUE, any numbers 
##' that are within eps of the default strategy are not replaced.
##' @export
K93_hyperpar <- make_K93_hyperpar()
