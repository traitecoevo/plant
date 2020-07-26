# Built from  R/ff16.R on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16
## We can probably actually do better than this with an S3 method on
## the actual strategy?  That would need to be organised by the
## templating though and that's stretched to the limit.

##' Create a K93 Plant or Cohort
##' @title Create a K93 Plant or Cohort
##' @param s A \code{\link{K93_Strategy}} object
##' @export
##' @rdname K93
##' @examples
##' pl <- K93_Plant()
##' pl$height
K93_Plant <- function(s=K93_Strategy()) {
  Plant("K93", "K93_Env")(s)
}

#' Compute the whole plant light compensation point for a single
#' plant with K93 strategy. Called via general function in plant.R
##' @export
##' @rdname K93
`lcp_whole_plant.Plant<K93>` <- function(p, ...) {
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
  Parameters("K93","K93_Env")(hyperpar=K93_hyperpar)
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

# ##' @export
# ##' @rdname K93
# K93_PlantPlus <- function(s=K93_Strategy()) {
#   PlantPlus("K93")(s)
# }

## Helper:
##' @export
##' @rdname K93_Environment
##' @param p A Parameters object
K93_make_environment <- function(p) {
  K93_Environment(p$disturbance_mean_interval, p$seed_rain, p$k_I, p$control)
}

##' Construct a fixed environment for FF16 strategy
##'
##' @param e=1.0 Value of environment 
##' @param p A Parameters object
##' @param height_max = 150.0 maximum possible height in environment
##' @rdname FF16_Environment
##'
##' @export
K93_fixed_environment <- function(e=1.0, p = K93_Parameters(), height_max = 300.0) {
  env <- K93_make_environment(p)
  env$set_fixed_environment(e, height_max)
  env
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

  function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }
    b_0 <- with_default("b_0")
    b_1 <- with_default("b_1")
    b_2 <- with_default("b_2")
    c_0 <- with_default("c_0")
    c_1 <- with_default("c_1")
    d_0 <- with_default("d_0")
    d_1 <- with_default("d_1")


    extra <- cbind(NULL)

    overlap <- intersect(colnames(m), colnames(extra))
    if (length(overlap) > 0L) {
      stop("Attempt to overwrite generated parameters: ",
           paste(overlap, collapse=", "))
    }

    ## Filter extra so that any column where all numbers are with eps
    ## of the default strategy are not replaced:
    if (filter) {
      if (nrow(extra) == 0L) {
        extra <- NULL
      } else {
        pos <- diff(apply(extra, 2, range)) == 0
        if (any(pos)) {
          eps <- sqrt(.Machine$double.eps)
          x1 <- extra[1, pos]
          x2 <- unlist(s[names(x1)])
          drop <- abs(x1 - x2) < eps & abs(1 - x1/x2) < eps
          if (any(drop)) {
            keep <- setdiff(colnames(extra), names(drop)[drop])
            extra <- extra[, keep, drop=FALSE]
          }
        }
      }
    }

    if (!is.null(extra)) {
      m <- cbind(m, extra)
    }
    m
  }

}

##' Hyperparameter function for K93 physiological model
##' @title Hyperparameter function for K93 physiological model
##' @export
K93_hyperpar <- make_K93_hyperpar()
