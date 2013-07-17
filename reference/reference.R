## As with reference.R, but for the simpler case of a single species.
source(system.file("tests/helper-tree.R", package="tree"))
source("reference-fun.R")

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.1)))
p$seed_rain <- 1.1
p$set_parameters(list(patch_area=1.0)) # See issue #13

path <- "ref-single"
## Evolve path hard coded to my machine for now.
## Running with version: d649534d6d92fd7343a57c9d26d51eec8e155f2d
path.evolve <- "~/Documents/Projects/veg/falster-traitdiversity/src"
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
              has_hash("055af115ee5e17fceb8dad0ea337e9a7a528aecb"))
})

ref <- tidyup.reference(output, p[[1]])

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

## Now, do a little post-processing.
res.1$values$mass.leaf <-
  translate(res.1$values$height, height.to.mass.leaf, p[[1]])
res.1$values$density <- exp(res.1$values$log.density)
res.1$values$density.mass.leaf <-
  res.1$values$density * dhdml(res.1$values$mass.leaf, p[[1]])
res.1$values$survival <- exp(-res.1$values$mortality)

## Good; these look about the same
matplot(res.1$time, res.1$values$height, type="l", col="black", lty=1,
        xlab="Time", ylab="Leaf mass")
matlines(ref$time,  ref$values$height,   type="l", col="grey", lty=1)

## Plot the top cohort:
plot(res.1$time, res.1$values$height[,1], type="l",
     xlab="Time", ylab="Leaf mass")
lines(ref$time,  ref$values$height[,1],   col="grey")

ref.sub <- resample.reference(ref, res.1$time)

## Difference and relative difference are actually quite good (about
## 1e-3 for relative difference).
plot(res.1$time,
     res.1$values$height[,1] - ref.sub$values$height[,1],
     xlab="Time", ylab="Difference in height")
plot(res.1$time,
     (res.1$values$height[,1] - ref.sub$values$height[,1]) /
     res.1$values$height[,1],
     xlab="Time", ylab="Relative difference in height")

## Expect agreement of the top cohort to .2%
expect_that(res.1$values$height[,1],
            equals(ref.sub$values$height[,1], tolerance=0.0002))

## Then, deal with the fact that we have a different set of data
## reporting, missing some values.  These should all be brand new
## cohorts.  It might be a good idea to rip these out during
## processing, or add them in in my version.
m <- ref.sub$values$height
m <- m[,-c(ncol(m), ncol(m)-1)]
i <- which(is.na(res.1$values$height) & !is.na(m))
expect_that(length(unique(m[i])), equals(1))
expect_that(unique(m[i]),
            equals(new(CohortTop, p[[1]])$height))
m[i] <- NA

## Top cohorts should agree strongly, the others less so.
expect_that(res.1$values$height[,1:10],
            equals(m[,1:10], tolerance=0.001))
expect_that(res.1$values$height, equals(m, tolerance=0.02))

## Mortalities look different, even for the top cohort.  So probably
## some work to do here.
matplot(res.1$time, res.1$values$mortality, type="l", col="black", lty=1,
        xlab="Time", ylab="Mortality")
matlines(ref$time,  ref$values$mortality,   type="l", col="grey", lty=1)

## here is the mortality variable for the top cohort.  We're starting
## at the same point, but growing at quite a different rate.
plot(res.1$time, res.1$values$mortality[,1], type="l")
lines(ref$time,  ref$values$mortality[,1],   col="grey")

## This makes me wonder if we're actually looking at different light
## environments [checked -- not the case]
plot(res.1$time,
     (res.1$values$mortality[,1] - ref.sub$values$mortality[,1]) /
     res.1$values$mortality[,1],
     xlab="Time", ylab="Relative difference in mortality")
plot(res.1$time,
     abs(res.1$values$mortality[,1] - ref.sub$values$mortality[,1]) /
     res.1$values$mortality[,1],
     xlab="Time", ylab="Relative difference in mortality", log="xy")

## Seed output:
## These look like they might be scaled incorrectly in my case.
matplot(res.1$time, res.1$values$seeds, type="l", col="black", lty=1,
        xlab="Time", ylab="Seeds")
matlines(ref$time,  ref$values$seeds,   type="l", col="grey", lty=1)

## ...but the difference is larger than that.  .75 looks possibly
## related to probability of survival during dispersal.
plot(res.1$time,
     (res.1$values$seeds[,1] - ref.sub$values$seeds[,1]) /
     res.1$values$seeds[,1],
     xlab="Time", ylab="Relative difference in seeds")

## Appling that fix brings this right down to about the right point,
## with a consistent bias, but low enough that it doesn't monsterously
## bother me.  The wavy line part is caused by a similar shape in the
## size plot.  The rest of the bias could be anything.
##
## I think I'm actually in favour of running this with the Pi_0 *not*
## calculated through in the CohortTop, because it can't possibly know
## that number.  We just need to make sure that we deal with it
## later.  So during post-processing make sure that calculation
## happens.
plot(res.1$time,
     (res.1$values$seeds[,1]*.25 - ref.sub$values$seeds[,1]) /
     res.1$values$seeds[,1]*.25,
     xlab="Time", ylab="Relative difference in seeds")

## Density:
matplot(res.1$time, res.1$values$log.density, type="l", col="black", lty=1,
        xlab="Time", ylab="Log density")
matlines(ref$time,  ref$values$log.density,   type="l", col="grey", lty=1)

matplot(res.1$time, res.1$values$density, type="l", col="black", lty=1,
        xlab="Time", ylab="Density")
matlines(ref$time,  ref$values$density,   type="l", col="grey", lty=1)

idx <- 50
plot(res.1$values$log.density[idx,])
plot(ref.sub$values$log.density[idx,])

## These look ok
plot(output$strategies[[1]]$bound_n[60,])
plot(log(output$strategies[[1]]$bound_n[60,]))
plot(log(ref$values$density.mass.leaf[60,]))
plot(log(ref.sub$values$density.mass.leaf[50,]))

plot(ref.sub$values$density.mass.leaf[50,1:50] /
     exp(res.1$values$density[50,1:50]))

plot(ref.sub$values$density.mass.leaf[50,1:50] /
     exp(res.1$values$density[50,1:50]), log="xy")
## This is not related to the difference in form
r <- unname(dmldh(res.1$values$height[50,1:50], p[[1]]))

## These actually all look really good.
idx <- 98
plot(ref.sub$light.env[[idx]])
lines(res.1$light.env[[idx]])
