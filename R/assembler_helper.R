##' Wrapper around assembler that tries to make things a bit easier.
##'
##' If the strategy is naive, then
##' @title Assembler Helper
##' @param sys0 Initial state of the system
##' @param births Strategy for births (naive or sampling)
##' @param ... Additional arguments passed to either
##' \code{\link{assembler_stochastic_naive}} or to
##' \code{\link{assembler_sample_positive}}.  In particular,
##' \code{compute_viable_fitness}, \code{run_type},
##' \code{jump_to_attractor} are supported by both and
##' \code{n_mutants}, \code{n_immigrants} and \code{n_sample} are
##' supported by the underlying birth functions.
##' @param vcv_p Mutantional variance as a fraction of the landscape
##' bounds
##' @author Rich FitzJohn
##' @export
assembler_helper <- function(sys0, births="naive", ..., vcv_p=0.001) {
  births <- match.arg(births, c("naive", "sampling"))
  if (births == "naive") {
    ## The mutational vcv needs to depend in the incoming bounds,
    ## rather than the ones we get after finding the viable space.
    vcv <- mutational_vcv_proportion(sys0$bounds, vcv_p)
    assembler_stochastic_naive(sys0, vcv, ...)
  } else if (births == "sampling") {
    assembler_sample_positive(sys0, ...)
  } else {
    stop("Unknown births type ", births)
  }
}
