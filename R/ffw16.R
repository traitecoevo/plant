## We can probably actually do better than this with an S3 method on
## the actual strategy?  That would need to be organised by the
## templating though and that's stretched to the limit.

##' Create a FFW16 Plant
##' @title Create a FFW16 Plant
##' @param s A \code{\link{FFW16_Strategy}} object
##' @export
##' @examples
##' pl <- FFW16_Plant()
##' pl$height
FFW16_Plant <- function(s=FFW16_Strategy()) {
  Plant("FFW16")(s)
}
