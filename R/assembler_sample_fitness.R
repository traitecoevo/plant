##' @export
assembler_sample_positive <- function(community0,
                                      n_sample=1L,
                                      seed_rain_eps=1e-3,
                                      compute_viable_fitness=TRUE,
                                      jump_to_attractor=FALSE,
                                      filename=NULL) {
  births_sys <- make_births_sample_positive(n_sample)
  deaths_sys <- make_deaths_stochastic_naive(seed_rain_eps)
  assembler(community0, births_sys, deaths_sys, filename,
            compute_viable_fitness, jump_to_attractor)
}

make_births_sample_positive <- function(n) {
  function(sys) {
    if (is.null(sys$bounds)) {
      ## No viable region
      ret <- matrix(nrow=0, ncol=1)
      colnames(ret) <- sys$trait_names
    } else {
      f <- fitness_landscape_approximate(sys)
      ret <- cbind(rejection_sample(n, f, sys$bounds))
      colnames(ret) <- sys$trait_names
      attr(ret, "landscape_approximate") <- f
    }
    ret
  }
}

## Quickly compute a spline function with n points based on the
## fitness function:
fitness_landscape_grid <- function(community, n=50,
                                   finite_only=TRUE, log_space=TRUE,
                                   force=TRUE, bounds=NULL) {
  if (!inherits(community, "community")) {
    stop("Expected a community object")
  }
  if (is.null(bounds)) {
    bounds <- community$bounds
  }
  if (nrow(bounds) != 1) {
    stop("Only working for one trait at the moment")
  }
  if (log_space) {
    x <- seq_log(bounds[[1]], bounds[[2]], n)
  } else {
    x <- seq(bounds[[1]], bounds[[2]], length.out=n)
  }
  mutant_seed_rain <- community$make_landscape(force)
  if (is.null(mutant_seed_rain)) {
    stop("Constructing fitness landscape failed")
  }
  m <- cbind(trait=x, fitness=log(mutant_seed_rain(x)))
  if (finite_only) {
    m <- m[is.finite(m[,2]),,drop=FALSE]
  }
  cbind(m)
}

##' @export
fitness_landscape_approximate <- function(community, n=50L,
                                          log_space=TRUE, force=TRUE,
                                          bounds=NULL) {
  xy <- fitness_landscape_grid(community, n, log_space=log_space,
                               force=force, bounds=bounds)
  splinefun_log(xy[,1], xy[,2])
}
