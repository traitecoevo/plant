assembler_sample_positive <- function(community0,
                                      n_sample=1L,
                                      seed_rain_eps=1e-3,
                                      compute_viable_fitness=FALSE,
                                      filename=NULL) {
  births_sys <- make_births_sample_positive(n_sample)
  ## TODO: Move into the core assembly algorithm.
  if (compute_viable_fitness) {
    message("Computing viable bounds")
    if (nrow(bounds) != 1) {
      stop("This is not going to work with multiple traits yet")
    }
    bounds <- viable_fitness(community0$trait_names,
                             community0$to_parameters(),
                             bounds=drop(bounds))
  }
  deaths_sys <- make_deaths_stochastic_naive(seed_rain_eps)
  assembler(community0, births_sys, deaths_sys, filename)
}

make_births_sample_positive <- function(n) {
  function(sys) {
    ## This will fail on the first step, because
    ## fitness_landscape_approximate does not work.
    f <- fitness_landscape_approximate(sys)
    rejection_sample(n, f, sys$bounds)
  }
}

## Quickly compute a spline function with n points based on the
## fitness function:
fitness_landscape_grid <- function(community, n=50,
                                   finite_only=TRUE, log_space=TRUE) {
  if (!inherits(community, "community")) {
    stop("Expected a community object")
  }
  bounds <- community$bounds
  if (log_space) {
    x <- seq_log(bounds[[1]], bounds[[2]], n)
  } else {
    x <- seq(bounds[[1]], bounds[[2]], length.out=n)
  }
  mutant_seed_rain <- community$make_landscape()
  m <- cbind(trait=x, fitness=log(mutant_seed_rain(x)))
  if (finite_only) {
    m <- m[is.finite(m[,2]),,drop=FALSE]
  }
  m
}

fitness_landscape_approximate <- function(sys, n=50L) {
  xy <- fitness_landscape_grid(sys, n)
  splinefun_log(xy[,1], xy[,2])
}
