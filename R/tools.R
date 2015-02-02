##' Create a list of Strategies or Plants by varying a single trait.
##'
##' @title Create a list of Strategies
##' @param x Values for the trait.  This must be a \emph{matrix}, with
##' column names corresponding to entries in \code{Strategy} and rows
##' representing different values.
##' @param strategy Default strategy to modify
##' @export
strategy_list <- function(x, strategy=Strategy()) {
  trait_names <- colnames(x)
  f <- function(xi) {
    strategy[trait_names] <- xi
    strategy
  }
  lapply(matrix_to_list(x), f)
}

##' @rdname strategy_list
##' @export
plant_list <- function(x, strategy=Strategy()) {
  lapply(strategy_list(x, strategy), Plant)
}

##' Helper function to create trait matrices suitable for
##' \code{\link{strategy_list}}.
##'
##' @title Create trait matrix
##' @param x Values
##' @param trait_name Name of a single trait
##' @export
##' @author Rich FitzJohn
trait_matrix <- function(x, trait_name) {
  m <- matrix(x, ncol=1)
  colnames(m) <- trait_name
  m
}
