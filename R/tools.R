##' Create a list of Strategies or Plants by varying a single trait.
##'
##' @title Create a list of Strategies
##' @param x Values for the trait
##' @param trait_name Name of the trait to vary (may be any element of
##' Strategy)
##' @param strategy Default strategy to modify
##' @export
strategy_list <- function(x, trait_name, strategy=Strategy()) {
  f <- function(xi) {
    strategy[trait_name] <- xi
    strategy
  }
  lapply(x, f)
}

##' @rdname strategy_list
plant_list <- function(x, trait_name, strategy=Strategy()) {
  lapply(strategy_list(x, trait_name, strategy), Plant)
}
