## We can probably actually do better than this with an S3 method on
## the actual strategy?  That would need to be organised by the
## templating though and that's stretched to the limit.

##' Create a FFW16 Plant or Cohort
##' @title Create a FFW16 Plant or Cohort
##' @param s A \code{\link{FFW16_Strategy}} object
##' @export
##' @rdname FFW16
##' @examples
##' pl <- FFW16_Plant()
##' pl$height
FFW16_Plant <- function(s=FFW16_Strategy()) {
  Plant("FFW16")(s)
}

##' @export
##' @rdname FFW16
FFW16_Cohort <- function(s=FFW16_Strategy()) {
  Cohort("FFW16")(s)
}

##' @export
##' @rdname FFW16
FFW16_Species <- function(s=FFW16_Strategy()) {
  Species("FFW16")(s)
}

##' @export
##' @rdname FFW16
##' @param ... Arguments!
FFW16_Parameters <- function(...) {
  Parameters("FFW16")(...)
}

##' @export
##' @rdname FFW16
##' @param p A \code{Parameters<FFW16>} object
FFW16_Patch <- function(p) {
  Patch("FFW16")(p)
}

##' @export
##' @rdname FFW16
FFW16_EBT <- function(p) {
  EBT("FFW16")(p)
}

##' @export
##' @rdname FFW16
FFW16_StochasticSpecies <- function(s=FFW16_Strategy()) {
  StochasticSpecies("FFW16")(s)
}

##' @export
##' @rdname FFW16
FFW16_StochasticPatch <- function(p) {
  StochasticPatch("FFW16")(p)
}

##' @export
##' @rdname FFW16
FFW16_StochasticPatchRunner <- function(p) {
  StochasticPatchRunner("FFW16")(p)
}

## Helper:
## TODO: consider directly using the C++ version in environment.h
##' @export
##' @rdname Environment
##' @param p A Parameters object
make_environment <- function(p) {
  Environment(p$disturbance_mean_interval,
              p$seed_rain,
              p$control)
}
