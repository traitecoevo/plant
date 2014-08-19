## This is a bit of a hack until we can get proper serialisation built
## into the objects.  I'm reluctant to do too much of that, because
## many of the objects (ebt and down) are not really designed to be
## serialised.

## Case getting the community load/save working.
library(tree)

ebt_base_parameters <- tree:::ebt_base_parameters

## From the successional diversity project:
assembler_lma_parameters <- function(time_disturbance, slope) {
  ## Specialised parameters:
  model <- list(c_r1=0.5,
                c_r2=0,
                B4=slope,
                ## keeps mean centered at global mean:
                a4=10^(0.1369 + slope*-0.9819))
  ## Build up a parameters object from all of that
  p <- ebt_base_parameters()
  p$strategy_default <- new(Strategy, model)
  p$disturbance <- new(Disturbance, time_disturbance)
  p
}

## Data members:
##   Contents of "$parameters"
##   control
##   [strategies] (will refuse to serialise if size > 0)
##   [seed_rain] (ditto)
##   [is_resident] (ditto)
##   disturbance (needs care)
##   strategy_default (needs care)
p <- assembler_lma_parameters(4, 1.2)
obj <- serialise_parameters(p)
p2 <- unserialise_parameters(obj)
obj2 <- serialise_parameters(p2)
all.equal(obj, obj2)
