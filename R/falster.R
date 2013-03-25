## My interface to Daniel's *R* version of the growth model.  This is
## used to cross check the C++ version.

## Daniel's version makes extensive use of global variables for both
## parameters and state.  So that I can more easily use this to test
## against my code without clobbering anything.  The basic idea will
## be to read things into an environment.

## The files falster/params.r and falster/growthModel.r are unchanged
## from versions recieved about 2012-10-29.
make.falster.generator <- function(path) {
  force(path)
  function() {
    e <- new.env()
    source(file.path(path, "falster/params.r"), local=e)
    source(file.path(path, "falster/growthModel.r"), local=e)
    ## Not used.
    ## source("falster/plots.r", local=e)
    ## TODO: These should be settable.
    e$traits <- list(lma=1.11E-01, rho=608, hmat=20, s=3.8e-5)

    ## Reset a few parameters to revert back to C++ defaults
    e$traits$lma <- 0.1978791
    e$traits$hmat <- 16.5958691
    e$p.c_bio <- 2.45e-2
    e$p.a4 <- 0.0286
    e$p.theta <- 4669

    ## Function to return current values of parameters.
    ## All parameters can be translated as x -> p.x, where 'x' is the
    ## name in the C++ version.
    vars <- c("eta", "theta", "a1", "B1", "a2", "B2", "a3", "b",
              "n_area", "c_Rl", "c_Rr", "c_Rs", "c_Rb", "Y",
              "c_bio", "a4", "B4", "k_b", "k_r", "c_p1", "c_p2",
              "c_acc", "c_d0", "c_d1", "c_d2", "c_d3", "c_r1", "c_r2")
    vars.msg <- list("Pi_0"=NA, "c_s0"=NA)
    get.pars <- function() {
      vars.traits <- list(lma=e$traits$lma, rho=e$traits$rho,
                          hmat=e$traits$hmat, s=e$traits$s)
      ret <- lapply(sprintf("p.%s", vars), get, e)
      names(ret) <- vars
      c(vars.traits, ret, vars.msg)
    }
    e$get_parameters <- get.pars

    ## These functions were not included in the core model, so I'm
    ## including them so that changes in the C++ version will be
    ## detected.
    source(file.path(path, "extra.R"), local=e)

    e
  }
}

make.falster <- make.falster.generator(getwd())
rm(make.falster.generator)
