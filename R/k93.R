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
    with_default <- function(name, default_value=s$pars[[name]]) {
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

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
K93_expand_state <- function(results) {

  # currently not doing anything
  
  results
}
