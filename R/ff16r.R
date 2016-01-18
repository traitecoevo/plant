## We can probably actually do better than this with an S3 method on
## the actual strategy?  That would need to be organised by the
## templating though and that's stretched to the limit.

##' Create a FF16r Plant or Cohort
##' @title Create a FF16r Plant or Cohort
##' @param s A \code{\link{FF16r_Strategy}} object
##' @export
##' @rdname FF16r
##' @examples
##' pl <- FF16r_Plant()
##' pl$height
FF16r_Plant <- function(s=FF16r_Strategy()) {
  Plant("FF16r")(s)
}

##' @export
##' @rdname FF16r
FF16r_Cohort <- function(s=FF16r_Strategy()) {
  Cohort("FF16r")(s)
}

##' @export
##' @rdname FF16r
FF16r_Species <- function(s=FF16r_Strategy()) {
  Species("FF16r")(s)
}

##' @export
##' @rdname FF16r
##' @param ... Arguments!
FF16r_Parameters <- function(...) {
  Parameters("FF16r")(...)
}

##' @export
##' @rdname FF16r
##' @param p A \code{Parameters<FF16r>} object
FF16r_Patch <- function(p) {
  Patch("FF16r")(p)
}

##' @export
##' @rdname FF16r
FF16r_SCM <- function(p) {
  SCM("FF16r")(p)
}

##' @export
##' @rdname FF16r
FF16r_StochasticSpecies <- function(s=FF16r_Strategy()) {
  StochasticSpecies("FF16r")(s)
}

##' @export
##' @rdname FF16r
FF16r_StochasticPatch <- function(p) {
  StochasticPatch("FF16r")(p)
}

##' @export
##' @rdname FF16r
FF16r_StochasticPatchRunner <- function(p) {
  StochasticPatchRunner("FF16r")(p)
}

##' @export
##' @rdname FF16r
FF16r_PlantPlus <- function(s=FF16r_Strategy()) {
  PlantPlus("FF16r")(s)
}

##' @title Hyperparameters for FF16r physiological model are same as FF16
##' @rdname FF16r_hyperpar
##' @export
FF16r_hyperpar <- make_FF16_hyperpar()
