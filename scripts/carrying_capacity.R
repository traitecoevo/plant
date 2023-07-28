
##' Compute the carrying capacity (equilibrium per-capita seed
##' production) for a set of values of a trait.  Each is considered in
##' isolation.
##'
##' @title Carrying Capacity
##' @param trait_matrix A matrix of traits
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param seed_rain Initial seed rain (optional)
##' @param parallel Use multiple processors?
##' @author Rich FitzJohn
##' @export
carrying_capacity <- function(trait_matrix, p, birth_rate=1,
                              parallel=FALSE) {
  f <- function(x) {
    ## NOTE: This is not very pretty because matrix_to_list, which we
    ## use to iterate over this, does not preserve array-ness or
    ## names.
    p <- expand_parameters(trait_matrix(x, traits),
                           remove_residents(p))
    p$birth_rate <- birth_rate
    equilibrium_birth_rate(p)$birth_rate
  }
  traits <- colnames(trait_matrix)
  unlist(loop(matrix_to_list(trait_matrix), f, parallel=parallel))
}

plant_log_console()
p <- scm_base_parameters()
ctrl <- equilibrium_quiet(Control())


## Find the viable range for lma for the default parameters:
r <- viable_fitness(bounds(lma=c(1e-5, 1e3)), p)
lma_detail <- seq_log_range(r, 201)

## Construct a fitness landscape across this region showing it is the
## right set of bounds:
lma <- trait_matrix(seq_log_range(r, length.out=51), "lma")
w <- fundamental_fitness(lma, p)

plot(w ~ lma, log="x", type="l", ylim=range(0, w), las=1,
     xlab="lma", ylab="Max growth rate")
abline(h=0, lty=2, col="red")
lines(lma_detail, splinefun_log(lma, w)(lma_detail),
      col="#0000ff55", lwd=5)

## This one is heaps slower, so easier to work with a smaller set of
## points.
lma <- trait_matrix(seq_log_range(r, length.out=31), "lma")
K <- carrying_capacity(lma, p, parallel=TRUE)

plot(K ~ lma, log="x", type="l", ylim=range(0, K), las=1,
     xlab="lma", ylab="Carrying capacity")
abline(h=0, lty=2, col="red")
lines(lma_detail, splinefun_log(lma, K)(lma_detail),
      col="#0000ff55", lwd=5)
