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
fitness_landscape_gp <- function(community,
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
  siefecor::interp1d(src, n_initial, n_final, n_each, cost=cost)
}

.R6_landscape_approximate <- R6::R6Class(
  "landscape_approximate",
  public=list(
    type=NULL,
    X=NULL,
    y=NULL,
    initialize=function(community, type, ...) {
      self$type <- match.arg(type, c("naive", "gp"))
      switch(self$type,
             naive=self$initialize_naive(community, ...),
             gp=self$initialize_gp(community, ...),
             stop("Not implemented")) # can't get here
    },
    initialize_naive=function(community, ...) {
      xy <- fitness_landscape_grid(community, ...)
      self$X <- xy[,1]
      self$y <- xy[,2]
      private$f <- splinefun_log(self$X, self$y)
    },
    initialize_gp=function(community, ...) {
      gp <- fitness_landscape_gp(community, ...)
      self$X <- gp$X
      self$y <- gp$y
      private$f <- landscape_approximate_gp_helper(gp)
      private$data <- gp$save()
    },
    predict=function(x) {
      private$f(x)
    },
    restore=function() {
      if (identical(self$type, "gp")) {
        gp <- siefecor:::gpreg_restore(private$data)
        private$f <- landscape_approximate_gp_helper(gp)
      }
    }
    ),
  private=list(
    f=NULL,
    data=NULL))
  
landscape_approximate <- function(...) {
  .R6_landscape_approximate$new(...)  
}

## This is here to try and force the enclosing scope to not be to
## aggressive.  This can probably be achived by using local({}) but I
## never *quite* know when that's right.
landscape_approximate_gp_helper <- function(gp) {
  function(x) {
    if (!is.matrix(x)) {
      x <- matrix(x, ncol=1)
    }
    gp$predict(log(x))
  }
}
