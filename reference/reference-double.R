source(system.file("tests/helper-tree.R", package="tree")) # has_hash

show.comparison <- function(res, ref, v, i=1) {
  if (interactive()) {
    ylim <- range(res$species[[i]][v,,], ref$species[[i]][v,,],
                  na.rm=TRUE, finite=TRUE)
    matplot(res$time, res$species[[i]][v,,], type="l", col="black", lty=1,
            xlab="Time", ylab=v, ylim=ylim)
    matlines(ref$time, ref$species[[i]][v,,], col="red", lty=2)
  }
}

show.comparisons <- function(res, ref.full, v) {
  n <- length(res$species)
  op <- par(mfrow=c(n, 1), mar=c(3.1, 4.1, .5, .5))
  on.exit(par(op))
  for (i in seq_len(n))
    show.comparison(res, ref.full, v, i)
}

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.0648406, hmat=26.3098)))
p$add_strategy(new(Strategy, list(lma=0.1977910, hmat=27.8790)))
p$seed_rain <- c(78, 26.6)
p$set_control_parameters(fast.control())
p$set_parameters(list(patch_area=1.0)) # See issue #13

path <- "ref-double"
path.evolve <- locate.evolve()

## Run falster-traitdiversity
run.reference(path, p, path.evolve=path.evolve, verbose=FALSE)

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
              has_hash("a5257c7490aaf017ff4ed39d19061c4490205bac"))
})

## This will always pass if the above test passes, of course.
test_that("Loaded two strategies worth of output", {
  expect_that(length(output$strategies), equals(2))
  expect_that(length(output$seed_rain_out), equals(2))
})

## Process the reference output to try and get everything into the
## same basic format as we expect.
##
## The reference output exists at a superset of times that our
## output exists at, because it is sampled at every ODE step (whereas
## we only collect values at cohort introductions).  So, subsample the
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

## In contrast with reference-single, the agreement here is not great:
expect_that(res$species[[1]]["height",,1],
            equals(ref.sub$species[[1]]["height",,1],
                   tolerance=0.002))
expect_that(res$species[[2]]["height",,1],
            equals(ref.sub$species[[2]]["height",,1],
                   tolerance=0.01))

## All variables look similar, qualitatively, but we've got a bit of a
## difference to deal with.
show.comparisons(res, ref.full, "height")
show.comparisons(res, ref.full, "log.mortality")
show.comparisons(res, ref.full, "log.density")
show.comparisons(res, ref.full, "seeds")
