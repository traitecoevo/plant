## This is a small helper function for calling tree from Python.  As
## we develop the ideas better, this will change.  The current verison
## allows computing of the fitness landscape for a set of parameters,
## similar to the test case we've been talking about.
##
## "time_disturbance" and "slope" are "p1" and "p2" from the
## discussion on the whiteboard.
##
## lma is a vector of species traits (the X values)
##
## seed_rain is a vector of species abundance (the Y values) and must
## be the same length as lma.
##
## equilibrium is a logical scalar, and switches between computing the
## equilibrium Y (default) or setting things up to compute
## non-equilibrium fitnesses.
##
## In the case where equilibrium=TRUE, the seed rain values are just
## starting points and the system will iterate until (approximate)
## equilibrium.
##' @export
for_python <- function(time_disturbance, slope, lma, seed_rain,
                       equilibrium=TRUE) {
  if (length(seed_rain) != length(lma)) {
    stop("seed_rain must be same length as lma")
  }

  p0 <- tree:::ebt_base_parameters()
  p0$strategy_default <- new(Strategy, list(c_r1=0.5, c_r2=0, B4=slope))
  p0$disturbance <- new(Disturbance, time_disturbance)

  ## Bounds to compute f over (in lma space)
  max_bounds <- rbind(lma=10^c(-4, 2))
  sys0 <- community(p0, "lma", bounds=max_bounds)
  ok <- sys0$set_viable_bounds(find_max_if_negative=TRUE)

  for (i in seq_along(lma)) {
    sys0$add_species(species(list(lma=lma[[i]]), seed_rain[[i]]))
  }

  if (equilibrium) {
    sys0$run_to_equilibrium()
  } else {
    sys0$run()
    sys0$set_seed_rain(seed_rain)
  }
  sys0
}
