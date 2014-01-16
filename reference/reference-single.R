## As with reference.R, but for the simpler case of a single species.
source(system.file("tests/helper-tree.R", package="tree"))
if (!evolve.is.installed())
  install.evolve()

## Simple plotting
show.comparison <- function(res, ref, v, i=1) {
  if (interactive()) {
    matplot(res$time, res$species[[i]][v,,], type="l", col="black", lty=1,
            xlab="Time", ylab=v)
    matlines(ref$time, ref$species[[i]][v,,], col="red", lty=2)
  }
}

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.1)))
p$seed_rain <- 1.1
p$set_control_parameters(fast.control())

path <- "ref-single"

## Run falster-traitdiversity
run.reference(path, p, verbose=FALSE)

## Load all the reference output:
output <- load.reference.output(path)
test_that("Output contains correct parameters", {
  expect_that(tree:::reference.compare.parameters(output$parameters, p),
              is_true())
})

## This needs updating every time that either the
## falster-traitdiversity is updated to give slightly different
## output, or when we process it differently with
## `load.reference.output()`.  It's mostly here to guard against the
## former, but most commonly fails because of the latter.
test_that("Output is as expected (strict test)", {
  expect_that(output[names(output) != "parameters"],
              has_hash("d832943c6809ac7417ae9ab4bf7b57d4eb9940b2"))
})

## The reference output exists at a superset of times that our output
## exists at, because it is sampled at every ODE step (whereas we only
## collect values at cohort introductions).  So, subsample the
## reference output to the times that our output corresponds to.
ref.full <- tidyup.reference.output(output)
ref.sub  <- resample.reference.output(ref.full)

## Run tree
p$set_parameters(list(patch_area=1.0)) # See issue #13
sched <- schedule.from.reference(ref.full)

res <- run.ebt.collect(p, sched)

############################################################################

#### 1: Disturbance regime.

## TODO: What to check about gamma?
age <- output$patch_age

test_that("Output disturbance calculations match", {
  ## tree version:
  d <- new(Disturbance, p$parameters[["mean_disturbance_interval"]])
  f.survival <- Vectorize(function(x) d$pr_survival(x))
  f.density  <- Vectorize(function(x) d$density(x))

  scal <- integrate(f.survival, 0, Inf)$value
  expect_that(f.density(age$age), equals(age$density))
  expect_that(f.survival(age$age) / scal, equals(age$density))
})

#### 2. Actual output

## a. Time is correct
test_that("Subsampled times agree",
          expect_that(res$time, equals(ref.sub$time)))

## b. Height vs time
show.comparison(res, ref.full, "height")

## The top few cohorts should really agree (0.02%, but the overall
## matrix of results differs by a little more [0.2%])
expect_that(res$species[[1]]["height",,1],
            equals(ref.sub$species[[1]]["height",,1],
                   tolerance=0.0002))

## This will work fairly well, but there are definitely differences:
expect_that(res$species[[1]]["height",,],
            equals(ref.sub$species[[1]]["height",,],
                   tolerance=0.003))

## c. Mortality vs time
show.comparison(res, ref.full, "log.mortality")
## Top cohorts should agree fairly strongly, the others less so.  We
## disagree here to about 3%.
expect_that(res$species[[1]]["log.mortality",,1:10],
            equals(ref.sub$species[[1]]["log.mortality",,1:10],
                   tolerance=0.005))
expect_that(res$species[[1]]["log.mortality",,],
            equals(ref.sub$species[[1]]["log.mortality",,],
                   tolerance=0.03))

## d. Seed output (weighted by survival)

## NOTE: I'm not sure I'm happy about this; we look a bit off here;
## the top cohort is a little overproductive, and the bottom cohorts
## are appear depressed.
show.comparison(res, ref.full, "seeds")
expect_that(res$species[[1]]["seeds",,],
            equals(ref.sub$species[[1]]["seeds",,],
                   tolerance=0.01))

## e. Density

## For density, we still are not sure how to translate between the two
## different units.  Because everything else agrees, it follows that
## these should be OK.
show.comparison(res, ref.full, "log.density")
expect_that(res$species[[1]]["log.density",,],
            equals(ref.sub$species[[1]]["log.density",,],
                   tolerance=0.03))
expect_that(exp(res$species[[1]]["log.density",,]),
            equals(exp(ref.sub$species[[1]]["log.density",,]),
                   tolerance=0.03))

## 3. Light environments
## These are slightly tricky in that they need rescaling too.
tr.light.env <- function(a, b)
  cbind(height=b[,"height"],
        canopy.openness=spline(a[,"height"], a[,"canopy.openness"],
          xout=b[,"height"])$y)
f <- function(i, res, ref)
  tr.light.env(res$light.env[[i]], ref$light.env[[i]])
res.light.env.stretched <- lapply(seq_along(res$light.env), f,
                                  res, ref.sub)

## The tolerance varies with length of time; we do get progressively
## worse as differences accumulate:
expect_that(res.light.env.stretched[1:10],
            equals(ref.sub$light.env[1:10]))
expect_that(res.light.env.stretched[30:40],
            equals(ref.sub$light.env[30:40], tolerance=1e-5))
expect_that(res.light.env.stretched,
            equals(ref.sub$light.env, tolerance=0.02))

## 4. Overall fitness
##
## These do differ a bit, but that's OK given how differently they are
## calculated.
expect_that(res$fitness,
            equals(unname(output$seed_rain_out), tolerance=0.01))
