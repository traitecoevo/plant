# Built from  R/ff16.R on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16

##' Create a K93 Individual or Node
##' @title Create a K93 Individual or Node
##' @param s A \code{\link{K93_Strategy}} object
##' @export
##' @rdname K93_Individual
##' @examples
##' pl <- K93_Individual()
##' pl$height
K93_Individual <- function(s=K93_Strategy()) {
  Individual("K93", "K93_Env")(s)
}

##' @export
##' @rdname FF16_Parameters
K93_Parameters <- function() {
  Parameters("K93","K93_Env")()
}

##' @export
##' @rdname FF16_make_environment
K93_make_environment <- function(light_availability_spline_tol = 1e-4, 
                                 light_availability_spline_nbase = 17,
                                 light_availability_spline_max_depth = 16, 
                                 light_availability_spline_rescale_usually = TRUE) {
  
  # for reasons unknown, we can't add arguments to the K93 constructor
  # as it causes the FF16 StochasticPatch tests to fail 🙃  opted to hard-code
  # these defaults into the K93_Environment
  
  e <- K93_Environment()
  
  e$light_availability <- ResourceSpline(light_availability_spline_tol, 
                     light_availability_spline_nbase, 
                     light_availability_spline_max_depth, 
                     light_availability_spline_rescale_usually)
  
  return(e)
}

##' @rdname FF16_fixed_environment
##' @export
K93_fixed_environment <- function(e=1.0, height_max = 300.0, ...) {
  env <- K93_make_environment(...)
  env$set_fixed_environment(e, height_max)
  env
}

##' @rdname FF16_test_environment
##' @examples
##' environment <- plant:::K93_test_environment(10)
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

  ret <- K93_make_environment()
  ret$light_availability$spline <- interpolator
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
##' @inheritParams FF16_hyperpar
##' @export
K93_hyperpar <- make_K93_hyperpar()
