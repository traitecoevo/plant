library(tree)

## Start by simulating just one patch with one species and seed it
## with one individual.

p <- new(Parameters)
p$add_strategy(new(Strategy))
p$set_parameters(list(patch_area=100.0))

## patch <- new(Patch, p)
patch <- new(PatchC, p)

run <- function(patch, t.max, n.max, reset=TRUE) {
  if (reset) {
    patch$reset()
    patch$add_seedlings(1)
  }

  tt <- numeric(0)
  hh <- list()

  while (patch$time < t.max && patch$n_individuals < n.max) {
    patch$step()
    tt <- c(tt, patch$time)
    hh <- c(hh, patch$height)
    cat(sprintf("%2.5f: %d / %d\n",
                patch$time, patch$ode_size/3, patch$n_individuals))
  }

  list(t=tt, n=sapply(hh, length), h=hh)
}

## This basically won't complete with 80,000 individuals using a plain
## Patch<Plant>; it grinds to a halt eventually.  Surprisingly, our
## memory usage is not that high, with max memory usage not too bad.
## Looks like we're more limited by time than space.

## step_deterministic is taking almost all of the time, and
## unsurprisingly this is through Solver::step() and eventually
## Patch<Plant>::derivs, and in turn Patch<Plant>::ode_values_set.

## Within this, compute_vars_phys takes about 53% of the time while
## construct_spline takes 45%.  The remaining 2% looks to be involved
## in actually setting values to the masses of plants (which will do
## things like update the various size traits of the individuals).

## Most of the time in compute_vars_phys is doing the actual
## integration, though a fraction of this is creating the integrator
## (which seems like something we should only do once, really).
## Another 2% is taken up by doing gsl_integration_workspace_free,
## which will be triggered by the destructor.  So we can save 4% by
## having this be a property of a plant.

## The spline construction time is basically all sunk into
## Plant::leaf_area_above()
set.seed(1)
res <- run(patch, 40, 150000)

plot(res$t, res$n, type="s")
plot(rep(res$t, res$n), unlist(res$h), log="y", pch=".",
     xlab="", ylab="Leaf mass")
