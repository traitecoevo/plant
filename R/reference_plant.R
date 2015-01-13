##' Create an instance of reference version of growth model
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
make_reference_plant <- function() {
  path <- system.file("reference_plant", package="tree2")

  e <- new.env()
  source(file.path(path, "R/params.r"), local=e)
  source(file.path(path, "R/growthModel.r"), local=e)

  e$traits <- list(lma=1.11E-01, rho=608, hmat=20, s=3.8e-5)
  ## Reset a few parameters to revert back to C++ defaults
  e$traits$lma <- 0.1978791
  e$traits$hmat <- 16.5958691
  e$p.c_bio <- 2.45e-2
  e$p.k_l0 <- 0.4565855
  e$p.theta <- 4669

  ## Function to return current values of parameters.
  ## All parameters can be translated as x -> p.x, where 'x' is the
  ## name in the C++ version.
  vars <- c("eta", "theta", "a1", "B1", "a3", "b",
            "n_area", "c_Rl", "c_Rr", "c_Rs", "c_Rb", "Y",
            "c_bio", "k_l0", "B4", "k_b", "k_s0", "B5", "k_r", "c_p1", "c_p2",
            "c_acc", "B7", "c_d0", "c_d1", "B6", "c_d2", "c_d3", "c_s0",
            "c_r1", "c_r2", "lma_0", "rho_0", "hmat_0", "s_0",  "n_area_0")
  get.pars <- function() {
    vars.traits <- list(lma=e$traits$lma, rho=e$traits$rho,
                        hmat=e$traits$hmat, s=e$traits$s)
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
