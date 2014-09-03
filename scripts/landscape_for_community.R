library(tree)

## This comes from our current set of biological explorations:
assembler_lma_parameters <- function(time_disturbance, slope) {
  model <- list(c_r1=0.5, c_r2=0,
                B4=slope)
  p <- tree:::ebt_base_parameters()
  p$strategy_default <- new(Strategy, model)
  p$disturbance <- new(Disturbance, time_disturbance)
  p
}

## Species traits that are *almost* at equilibrium
traits <- cbind(lma=c(0.0123404040290072, 0.00281884510501954))
## Seed rain (~= species density)
seed_rain <- c(77.13723,287.0169)
## Base parameter set:
p0 <- assembler_lma_parameters(32, 1.05)

## Build a community containing these species
obj <- make_community(traits, seed_rain, p0)

## And run this to equilibrium seed rain (alternatively, obj$run()
## will just set things up but hold the seed rain the same as what was
## given).
obj$run_to_equilibrium()

## This makes the fitness function; it takes a vector of lma values
## and will return fitness.  Below zero: can't invade, above zero: can
## invade.  If you want to compute a bunch of values at once, pass
## them in at once and this goes a bit faster.
w <- obj$make_landscape()

## Fitness for 50 mutants spaced equally in log-lma space.  It takes a
## while to run, sorry.
x <- seq_log(0.001, 0.1, 50)
y <- w(x)

plot(x, y, log="x")
abline(h=0, lty=2)
