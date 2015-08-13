## We can probably actually do better than this with an S3 method on
## the actual strategy?  That would need to be organised by the
## templating though and that's stretched to the limit.

##' Create a FFdev Plant or Cohort
##' @title Create a FFdev Plant or Cohort
##' @param s A \code{\link{FFdev_Strategy}} object
##' @export
##' @rdname FFdev
##' @examples
##' pl <- FFdev_Plant()
##' pl$height
FFdev_Plant <- function(s=FFdev_Strategy()) {
  Plant("FFdev")(s)
}

##' @export
##' @rdname FFdev
FFdev_Cohort <- function(s=FFdev_Strategy()) {
  Cohort("FFdev")(s)
}

##' @export
##' @rdname FFdev
FFdev_Species <- function(s=FFdev_Strategy()) {
  Species("FFdev")(s)
}

##' @export
##' @rdname FFdev
##' @param ... Arguments!
FFdev_Parameters <- function(...) {
  Parameters("FFdev")(...)
}

##' @export
##' @rdname FFdev
##' @param p A \code{Parameters<FFdev>} object
FFdev_Patch <- function(p) {
  Patch("FFdev")(p)
}

##' @export
##' @rdname FFdev
FFdev_EBT <- function(p) {
  EBT("FFdev")(p)
}

##' @export
##' @rdname FFdev
FFdev_StochasticSpecies <- function(s=FFdev_Strategy()) {
  StochasticSpecies("FFdev")(s)
}

##' @export
##' @rdname FFdev
FFdev_StochasticPatch <- function(p) {
  StochasticPatch("FFdev")(p)
}

##' @export
##' @rdname FFdev
FFdev_StochasticPatchRunner <- function(p) {
  StochasticPatchRunner("FFdev")(p)
}

##' @export
##' @rdname FFdev
FFdev_PlantPlus <- function(s=FFdev_Strategy()) {
  PlantPlus("FFdev")(s)
}
