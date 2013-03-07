library(methods)
library(Rcpp)
library(testthat)


if ( packageVersion("Rcpp") < package_version("0.10.0") )
  error("Needs Rcpp at least version 0.10.0")

if ( !exists("Spline") ) {
  ## TODO: Sort out how we run this.  Package may help.
  dyn.load(file.path("../src", "tree.so"))
  tree_module <- Module("tree", "tree")
  
  Spline            <- tree_module$Spline
  AdaptiveSpline    <- tree_module$AdaptiveSpline
  RAdaptiveSpline   <- tree_module$RAdaptiveSpline
  
  .Call("set_sane_gsl_error_handling")
}
