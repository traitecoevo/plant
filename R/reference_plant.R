##' Create an instance of reference version for specified growth model
##'
##' @title Reference Growth Model
##' @return An object with a number of methods (accessed with
##' \code{$}).  These will hopefully be documented at some point.
##' @author Daniel S. Falster (extra code and porting by Rich
##' FitzJohn)
##' @param type name of model to make a reference plant for (only FF16
##' supported at present).
##' @export
make_reference_plant <- function(type="FF16") {
  if (type == "FF16") {
    ret <- make_reference_plant_FF16()
  } else {
    stop("No reference plant for type ", type)
  }
  ret
}

##' Create an instance of reference version of growth model FF16
##'
##' This is an old version of the growth model, implemented by Daniel
##' Falster.  It does \emph{not} implement any of the actual growth.
##' It is used in testing code to make sure that the full model
##' version gives sensible results.
##'
##' @title Reference Growth Model
##' @return An object with a number of methods (accessed with
##' \code{$}).  These will hopefully be documented at some point.
##' @author Daniel S. Falster (extra code and porting by Rich
##' FitzJohn)
##' @export
make_reference_plant_FF16 <- function() {

  path <- system.file("reference_plant_ff16", package=.packageName)

  e <- new.env()
  source(file.path(path, "R/params.r"), local=e)
  source(file.path(path, "R/growthModel.r"), local=e)

  e$traits <- list(lma=1.11E-01, rho=608, hmat=20, omega=3.8e-5)
  ## Reset a few parameters to revert back to C++ defaults
  e$traits$lma <- 0.1978791
  e$traits$hmat <- 16.5958691
  e$p.a_bio <- 2.45e-2
  e$p.k_l <- 0.4565855
  e$p.theta <- 4669

  ## Function to return current values of parameters.
  ## All parameters can be translated as x -> p.x, where 'x' is the
  ## name in the C++ version.
  vars <- c("eta", "theta", "a_l1", "a_l2", "a_r1", "a_b1",
            "n_area", "r_l", "r_r", "r_s", "r_b", "a_y",
            "a_bio", "k_l", "k_b", "k_s", "k_r", "a_p1", "a_p2",
            "a_f3","d_I", "a_dG1", "a_dG2", "a_d0",
            "a_f1", "a_f2")
  get.pars <- function() {
    vars.traits <- list(lma=e$traits$lma, rho=e$traits$rho,
                        hmat=e$traits$hmat, omega=e$traits$omega)
    ret <- lapply(sprintf("p.%s", vars), get, e)
    names(ret) <- vars
    c(vars.traits, ret)
  }
  e$get_parameters <- get.pars

  ## These functions were not included in the core model, so I'm
  ## including them so that changes in the C++ version will be
  ## detected.
  source(file.path(path, "extra.R"), local=e)

  e
}
