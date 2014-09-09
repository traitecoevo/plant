##' @export
assembler_sample_positive <- function(community0, n_sample=1L,
                                      seed_rain_eps=1e-3,
                                      approximate_type="naive",
                                      ...) {
  births_sys <- make_births_sample_positive(n_sample, approximate_type)
  deaths_sys <- make_deaths_stochastic_naive(seed_rain_eps)
  assembler(community0, births_sys, deaths_sys, ...)
}

make_births_sample_positive <- function(n, approximate_type="naive") {
  approximate_type <- match.arg(approximate_type, c("naive", "gp"))
  function(sys) {
    if (is.null(sys$bounds)) {
      ## No viable region
      ret <- matrix(nrow=0, ncol=1)
      colnames(ret) <- sys$trait_names
    } else {
      f <- fitness_landscape_approximate(sys, approximate_type)
      ret <- cbind(rejection_sample(n, f, sys$bounds))
      colnames(ret) <- sys$trait_names
      attr(ret, "landscape_approximate") <- f
      attr(ret, "done") <- nrow(ret) == 0
    }
    ret
  }
}

## Quickly compute a spline function with n points based on the
## fitness function:
fitness_landscape_grid <- function(community, n=50,
                                   finite_only=TRUE, log_space=TRUE,
                                   bounds=NULL) {
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
  mutant_fitness <- community$make_landscape()
  m <- cbind(trait=x, fitness=mutant_fitness(x))
  if (finite_only) {
    m <- m[is.finite(m[,2]),,drop=FALSE]
  }
  cbind(m)
}

## This should basically be enough.  Need to work out how to tidy up
## all the options, and to return a function if we're replacing the
## spline version.  Then work out how to nicely switch between the
## different approaches.

## This should basically be enough.  Need to work out how to tidy up
## all the options, and to return a function if we're replacing the
## spline version.  Then work out how to nicely switch between the
## different approaches.
fitness_landscape_approximate_gp <- function(community,
                                             n_initial=20, n_final=50,
                                             n_each=5, n_predict=500,
                                             log_space=TRUE,
                                             bounds=NULL,
                                             lower_limit=NULL) {
  loadNamespace("siefecor") # makes this optional
  if (!inherits(community, "community")) {
    stop("Expected a community object")
  }
  if (is.null(bounds)) {
    bounds <- community$bounds
  }
  if (nrow(bounds) != 1) {
    stop("Only working for one trait at the moment")
  }

  mutant_fitness <- community$make_landscape()
  if (log_space) {
    objective <- function(x) {
      mutant_fitness(exp(x))
    }
    bounds <- log(bounds)
  } else {
    objective <- mutant_fitness
  }
  cost <- siefecor::make_cost_var_scaled_capped(scal=10, mu_max=0.0)
  x <- cbind(seq(bounds[[1]], bounds[[2]], length.out=n_predict))
  src <- siefecor::data_source$new(x, objective,
                                   lower_limit=lower_limit,
                                   verbose=TRUE)
  res <- siefecor::interp1d(src, n_initial, n_final, n_each, cost=cost)
  function(x) {
    if (!is.matrix(x)) {
      x <- matrix(x, ncol=1)
    }
    res$predict(log(x))
  }
}

##' @export
fitness_landscape_approximate_naive <- function(community, n=50L,
                                          log_space=TRUE,
                                          bounds=NULL) {
  xy <- fitness_landscape_grid(community, n, log_space=log_space,
                               bounds=bounds)
  splinefun_log(xy[,1], xy[,2])
}

##' @export
fitness_landscape_approximate <- function(community, type="naive",
                                          ...) {
  switch(match.arg(type, c("naive", "gp")),
         naive=fitness_landscape_approximate_naive(community, ...),
         gp=fitness_landscape_approximate_gp(community, ...),
         stop("Can't get here"))
}
