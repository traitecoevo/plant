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
