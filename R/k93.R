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

##' Hyperparameters for K93 physiological model
##' @title Hyperparameters for K93 physiological model
##' @export
##' @rdname K93_hyperpar
make_K93_hyperpar <- function() {
}

##' Hyperparameter function for K93 physiological model
##' @title Hyperparameter function for K93 physiological model
##' @export
K93_hyperpar <- make_K93_hyperpar()
