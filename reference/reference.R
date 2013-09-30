## As with reference.R, but for the simpler case of a single species.
source(system.file("tests/helper-tree.R", package="tree"))
source("reference-fun.R")

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.1)))
p$seed_rain <- 1.1
p$set_parameters(list(patch_area=1.0)) # See issue #13

path <- "ref-single"
path.evolve <- locate.evolve()

## Run falster-traitdiversity
run.reference(path, p, path.evolve=path.evolve, verbose=FALSE)

## Load all the reference output:
output <- load.reference.output(path)
test_that("Output contains correct parameters", {
  expect_that(compare.parameters(output$parameters, p), is_true())
})

## This needs updating every time that either the
## falster-traitdiversity is updated to give slightly different
## output, or when we process it differently with
## `load.reference.output()`.  It's mostly here to guard against the
## former, but most commonly fails because of the latter.
test_that("Output is as expected (harsh test)", {
  expect_that(output[names(output) != "parameters"],
              has_hash("f6d22a6383273d99491feacb78cee2eb113f4186"))
})

## The reference output exists at a superset of times that our output
## exists at, because it is sampled at every ODE step (whereas we only
## collect values at cohort introductions).  So, subsample the
## reference output to the times that our output corresponds to.
ref.full <- tidyup.reference(output, p[[1]])
ref.sub  <- resample.reference(ref.full)

## 1: Disturbance:
## TODO: What to check about gamma?
age <- output$patch_age

test_that("Output disturbance calculations match", {
  ## tree version:
  d <- new(Disturbance, p$parameters[["mean_disturbance_interval"]])
  f.survival <- Vectorize(function(x) d$survival_probability(0, x))
  f.density  <- Vectorize(function(x) d$density(x))

  scal <- integrate(f.survival, 0, Inf)$value
  expect_that(f.density(age$age), equals(age$density))
  expect_that(f.survival(age$age) / scal, equals(age$density))
})

## Now, run tree with same parameters, but relaxed control parameters:
ctrl.new <- list()
ctrl.new$environment_light_rescale_usually <- TRUE
ctrl.new$environment_light_tol <- 1e-4
ctrl.new$plant_assimilation_rule <- integrator_gsl_rule("GAUSS21")
ctrl.new$plant_assimilation_over_distribution <- FALSE
ctrl.new$plant_assimilation_tol <- 1e-4
ctrl.new$ode_tol_rel <- 1e-4
ctrl.new$ode_tol_abs <- 1e-4
ctrl.new$ode_step_size_max <- 5
p$set_control_parameters(ctrl.new)

## Build a cohort schedule from the EBT's schedule:
sched <- new(CohortSchedule, p$size)
sched$set_times(output$patch_age$age[output$patch_age$add.cohort.1],
                1)
sched$max_time <- max(output$patch_age$age)
## Force ODE times to be the same as Daniel's
##
## NOTE: it would be nice to be able to *step* the system here, rather
## than doing advance.  And have it use the current set of ode times.
sched$ode_times <- output$patch_age$age

ebt <- new(EBT, p)
ebt$cohort_schedule <- sched

## Do the run
res.1 <- run.ebt(ebt)

## Now, do a little post-processing, adding on variables that we want
## to compare with the reference output.
res.1$values$mass.leaf <-
  translate(res.1$values$height, height.to.mass.leaf, p[[1]])
res.1$values$density <- exp(res.1$values$log.density)
res.1$values$density.mass.leaf <-
  res.1$values$density * dhdml(res.1$values$mass.leaf, p[[1]])
res.1$values$survival <- exp(-res.1$values$mortality)
## These are scaled differently:
##
## NOTE: I think I'm actually in favour of running this with the Pi_0
## *not* calculated through in the CohortTop, because it can't
## possibly know that number.  We just need to make sure that we deal
## with it later.  So during post-processing make sure that
## calculation happens.
res.1$values$seeds <- res.1$values$seeds * p$parameters$Pi_0

## Start the comparisons:

show.comparison <- function(res, ref.full, ref.sub, v) {
  if (interactive()) {
    matplot(res$time, res[[v]], type="l", col="black", lty=1,
            xlab="Time", ylab=paste(v, collapse=": "))
    matlines(ref.full$time, ref.full[[v]], type="l", col="red", lty=2)
    matlines(ref.sub$time, ref.sub[[v]], type="l", col="blue", lty=2)
  }
}

## a. Height vs time
show.comparison(res.1, ref.full, ref.sub, c("values", "height"))

## The top few cohorts should really agree (0.02%, but the overall
## matrix of results differs by a little more [0.2%])
expect_that(res.1$values$height[,1],
            equals(ref.sub$values$height[,1], tolerance=0.0002))
expect_that(res.1$values$height,
            equals(ref.sub$values$height, tolerance=0.002))

## b. Mortality
show.comparison(res.1, ref.full, ref.sub, c("values", "mortality"))

## Top cohorts should agree fairly strongly, the others less so.  We
## disagree here to about 3%.
expect_that(res.1$values$mortality[,1:10],
            equals(ref.sub$values$mortality[,1:10], tolerance=0.005))
expect_that(res.1$values$mortality,
            equals(ref.sub$values$mortality, tolerance=0.03))

## When translated into survival (which is what actually matters for
## the model) the disagreement is much better.
show.comparison(res.1, ref.full, ref.sub, c("values", "survival"))
expect_that(res.1$values$survival[,1:10],
            equals(ref.sub$values$survival[,1:10], tolerance=0.001))
expect_that(res.1$values$survival,
            equals(ref.sub$values$survival, tolerance=0.005))

## c. Seed output (weighted by survival)

## NOTE: I'm not sure I'm happy about this; we look a bit off here;
## the top cohort is a little overproductive, and the bottom cohorts
## are appear depressed.
show.comparison(res.1, ref.full, ref.sub, c("values", "seeds"))

## Agree to within 1%
expect_that(res.1$values$seeds,
            equals(ref.sub$values$seeds, tolerance=0.01))

## d. Density

## For density, we still are not sure how to translate between the two
## different units.  Because everything else agrees, it follows that
## these should be OK.

## 3. Light environments

## These are slightly tricky in that they need rescaling too.
tr.light.env <- function(a, b)
  cbind(height=b[,"height"],
        canopy.openness=spline(a[,"height"], a[,"canopy.openness"],
          xout=b[,"height"])$y)

f <- function(i)
  tr.light.env(res.1$light.env[[i]], ref.sub$light.env[[i]])
res.1.light.env.stretched <- lapply(seq_along(res.1$light.env), f)

## The tolerance varies with length of time; we do get progressively
## worse as differences accumulate:
expect_that(res.1.light.env.stretched[1:10],
            equals(ref.sub$light.env[1:10]))
expect_that(res.1.light.env.stretched[30:40],
            equals(ref.sub$light.env[30:40], tolerance=1e-5))
expect_that(res.1.light.env.stretched,
            equals(ref.sub$light.env, tolerance=0.02))

## 4. Overall fitness
fitness.ref <- unname(output$seed_rain_out)

## We need both cohort introduction times, and patch weights.  At the
## moment, this is a bit of a hack.
a <- sched$times(1)

## Here are the patch weights, for patches born at times 'a'.
d <- new(Disturbance, p$parameters[["mean_disturbance_interval"]])
pa <- sapply(a, function(ai) d$density(ai))

## Here is the total seed production -- this is per capita so we
## multiply by the incoming seed rain.
seeds <- res.1$values$seeds[nrow(res.1$values$seeds),]
fitness.res <- trapezium(a, pa * seeds) * p$seed_rain

## These do differ a bit, but that's OK given how differently they are
## calculated.
expect_that(fitness.res, equals(fitness.ref, tolerance=0.01))
