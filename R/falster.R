## My interface to Daniel's *R* version of the growth model.  This is
## used to cross check the C++ version.

## Daniel's version makes extensive use of global variables for both
## parameters and state.  So that I can more easily use this to test
## against my code without clobbering anything.  The basic idea will
## be to read things into an environment.

## The files falster/params.r and falster/growthModel.r are unchanged
## from versions recieved about 2012-10-29.
make.falster <- function() {
  e <- new.env()
  source("falster/params.r", local=e)
  source("falster/growthModel.r", local=e)
  ## Not used.
  ## source("falster/plots.r", local=e)
  e$traits <- list(lma=1.11E-01, rho=608, hmat=20)
  ## These are done in main.R, but skipping for now.
  ## e$p.c_r1 = 0.75
  ## e$p.c_r2 = 10

  ## The next thing that I need is a translation table for the different
  ## parameters.

  ## My version:
  e$tr.pars <-
    c("A0"=NA, # missing -- only used in LRC calc?
      "Y"="p.Y",
      "a1"="p.a1",
      "a2"="p.a2",
      "a3"="p.a3",
      "a4"="p.a4",
      "b"="p.b",
      "b1"="p.B1",
      "b2"="p.B2",
      "b4"="p.B4",
      "c_Rl"="p.c_Rl",
      "c_Rr"="p.c_Rr",
      "c_Rs"="p.c_Rs",
      "c_acc"="p.c_acc",
      "c_bio"="p.c_bio",
      "c_p1"="p.c_p1",
      "c_p2"="p.c_p2",
      "c_r1"="p.c_r1",
      "c_r2"="p.c_r2",
      "hmat"=NA,  # traits$hmat
      "kb"="p.k_b",
      "kr"="p.k_r",
      "lma"=NA, # traits$lma
      "ml0"=NA, # computed, from s
      "n"="p.eta",
      "nc"=NA, # I compute this -- not sure if Daniel does this too
      "rho"=NA, # traits$rho
      "s"=NA,   # traits$s
      "theta"="p.theta",
      "v"="p.n_area")
  e
}

## These (a5, B5) are used for computing the diameter (via Diameter).
## Not sure if/how these are used though.
## "p.a5" # not sure
## "p.B5" # not sure

## Other queries/comments
## c_r1, c_r2 (used in seed production) are altered in main.R to
## different values than we have.
##
## hmat appears twice in Daniels code, but once in mine.
