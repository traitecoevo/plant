# Built from  R/ff16.R on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16

## We can probably actually do better than this with an S3 method on
## the actual strategy?  That would need to be organised by the
## templating though and that's stretched to the limit.

##' Create a K93 Individual or Node
##' @title Create a K93 Individual or Node
##' @param s A \code{\link{K93_Strategy}} object
##' @export
##' @rdname K93
##' @examples
##' pl <- K93_Individual()
##' pl$height
K93_Individual <- function(s=K93_Strategy()) {
  Individual("K93", "K93_Env")(s)
}

##' @export
##' @rdname K93
K93_Node <- function(s=K93_Strategy()) {
  Node("K93", "K93_Env")(s)
}

##' @export
##' @rdname K93
K93_Species <- function(s=K93_Strategy()) {
  Species("K93", "K93_Env")(s)
}

##' @export
##' @rdname K93
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


## Helper to create K93_environment object. Useful for running individuals
##' create K93_environment object
##' @param shading_spline_tol
##'
##' @param shading_spline_nbase
##' @param shading_spline_max_depth
##' @param shading_spline_rescale_usually
##'
##' @export
##' @rdname K93_make_environment
K93_make_environment <- function(shading_spline_tol = 1e-4, 
                                 shading_spline_nbase = 17,
                                 shading_spline_max_depth = 16, 
                                 shading_spline_rescale_usually = TRUE) {
  
  # for reasons unknown, we can't add arguments to the K93 constructor
  # as it causes the FF16 StochasticPatch tests to fail 🙃  opted to hard-code
  # these defaults into the K93_Environment
  
  e <- K93_Environment()
  
  e$shading <- Shading(shading_spline_tol, 
                     shading_spline_nbase, 
                     shading_spline_max_depth, 
                     shading_spline_rescale_usually)
  
  return(e)
}

##' Construct a fixed environment for K93 strategy
##'
##' @param e Value of environment (default=1.0)
##' @param height_max = 300.0 maximum possible height in environment
##' @rdname K93_Environment
##'
##' @export
K93_fixed_environment <- function(e=1.0, height_max = 300.0) {
  env <- K93_make_environment()
  env$set_fixed_environment(e, height_max)
  env
}


##' This makes a pretend light environment over the plant height,
##' slightly concave up, whatever.
##' @title Create a test environment for K93 startegy
##' @param height top height of environment object
##' @param n number of points
##' @param light_env function for light environment in test object
##' @param n_strategies number of strategies for test environment
##' @export
##' @rdname K93_test_environment
##' @examples
##' environment <- K93_test_environment(10)
K93_test_environment <- function(height, n=101, light_env=NULL,
                             n_strategies=1) {
  hh <- seq(0, height, length.out=n)
  if (is.null(light_env)) {
    light_env <- function(x) {
      # arbitary function. aiming to produce values of -log(light_env)/0.01
      # in range 0:100
      exp(x/(height*2)) - (exp(.5) - 1)
    }
  }
  ee <- light_env(hh)
  interpolator <- Interpolator()
  interpolator$init(hh, ee)

  # parameters <- K93_Parameters()
  # parameters$strategies <- rep(list(K93_Strategy()), n_strategies)
  # 

  ret <- K93_make_environment()
  ret$shading$spline <- interpolator
  attr(ret, "light_env") <- light_env
  ret
}

##' Construct hyperparameter object for K93 physiological model
##' @title Hyperparameters for K93 physiological model
##' @param b_0 Growth intercept year-1
##' @param b_1 Growth asymptote year-1.(ln cm)-1
##' @param b_2 Growth suppression rate m2.cm-2.year-1
##' @param c_0 Mortality intercept year-1
##' @param c_1 Mortality suppression rate m2.cm-2.year-1
##' @param d_0 Recruitment rate (cm2.year-1)
##' @param d_1 Recruitment suppression rate (m2.cm-2)
##' @param eta Crown shape parameter
##' @param k_I Extinction coefficient used when estimating competitive effect
##' @export
##' @rdname make_K93_hyperpar
make_K93_hyperpar <- function(
        b_0 = 0.059,    # Growth intercept year-1
        b_1 = 0.012,    # Growth asymptote year-1.(ln cm)-1
        b_2 = 0.00041,  # Growth suppression rate m2.cm-2.year-1
        c_0 = 0.008,    # Mortality intercept year-1
        c_1 = 0.00044,  # Mortality suppression rate m2.cm-2.year-1
        d_0 = 0.00073,  # Recruitment rate (cm2.year-1)
        d_1 = 0.044,    # Recruitment suppression rate (m2.cm-2)
        eta = 12,       # Canopy shape parameter
        k_I = 0.01      # Scaling factor for competition
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
  assert_scalar(eta)
  assert_scalar(k_I)

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
##' @rdname K93_hyperpar
##' @export
K93_hyperpar <- make_K93_hyperpar()
