source(system.file("tests/helper-tree.R", package="tree"))

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.1)))
p$add_strategy(new(Strategy, list(lma=0.3)))
p$seed_rain <- c(1.1, 2.2)

path <- "ref"
## Evolve path hard coded to my machine for now.
## Running with version: d649534d6d92fd7343a57c9d26d51eec8e155f2d
path.evolve <- "~/Documents/Projects/veg/falster-traitdiversity/src"
## run.reference(path, p, path.evolve=path.evolve, verbose=TRUE)

output <- load.reference.output(path)
test_that("Output contains correct parameters", {
  compare.parameters <- function(p1, p2)
    isTRUE(all.equal(tree:::reference.from.parameters(p1),
                     tree:::reference.from.parameters(p2)))
  expect_that(compare.parameters(output$parameters, p), is_true())
})

## This needs updating every time that either the falster-growthmodel
## is updated to give slightly different output, or when we process it
## differently with `load.reference.output()`.  It's mostly here to
## guard against the former, but most commonly fails because of the
## latter.
test_that("Output is as expected (harsh test)", {
  expect_that(output[names(output) != "parameters"],
              has_hash("924331f149ba4e921a9f5bc79534a439be4a089c"))
})

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

## 2: Find out when the cohorts were introduced.

## Build a schedule.
sched <- new(CohortSchedule, p$size)
sched$set_times(output$patch_age$age[output$patch_age$add.cohort.1],
                1)
sched$set_times(output$patch_age$age[output$patch_age$add.cohort.2],
                2)
sched$max_time <- max(output$patch_age$age)

## Relaxed control parameters
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

ebt <- new(EBT, p)
ebt$cohort_schedule <- sched

## This seems like a function we're trying to use a lot.  Push into a
## general library of functions?
run.ebt <- function(ebt) {
  ## See test-ebt.R
  list.to.matrix <- function(x) {
    n <- max(sapply(x, length))
    t(sapply(x, function(i) c(i, rep(NA, n-length(i)))))
  }

  tt <- hh <- NULL
  ebt$reset()

  while (ebt$cohort_schedule$remaining > 0) {
    ebt$run_next()
    tt <- c(tt, ebt$time)
    hh <- c(hh, list(ebt$patch$height))
  }

  hh <- apply(matrix(unlist(hh, FALSE), p$size),
              1, list.to.matrix)
  list(time=tt, height=hh)
}

## Takes 1.4s
res.1 <- run.ebt(ebt)

## Now, set the ode times to match Daniel's
sched$ode_times <- output$patch_age$age
ebt$reset()
ebt$cohort_schedule <- sched

## Takes 1.2s
res.2 <- run.ebt(ebt)

expect_that(res.2, equals(res.1, tolerance=1e-5))

matplot(res.1$time, res.1$height[[1]], type="l", col="black", lty=1,
        xlab="Time", ylab="Height")

## This all suggests that we need to be quicker.  This translation
## takes 0.5s!  Not sure what the nicest solution is, but this is
## clearly not sustainable.
height.to.mass.leaf <- function(height, strategy) {
  p <- new(Plant, strategy)
  f <- function(x) {
    if (!is.finite(x)) return(NA)
    p$height <- x
    p$vars_size[["mass_leaf"]]
  }
  sapply(height, f)
}

mass.leaf.to.height <- function(mass.leaf, strategy) {
  pars <- strategy$parameters
  a1 <- pars$a1
  B1 <- pars$B1
  lma <- pars$lma
  a1 * pow(mass.leaf / lma, B1)
}
height.to.mass.leaf <- function(height, strategy) {
  pars <- strategy$parameters
  a1 <- pars$a1
  B1 <- pars$B1
  lma <- pars$lma
  (height / a1) ^ (1 / B1) * lma
}

res.1$mass.leaf <- lapply(res.1$height, function(x)
                          array(height.to.mass.leaf(x, p[[1]]), dim(x)))
res.2$mass.leaf <- lapply(res.2$height, function(x)
                          array(height.to.mass.leaf(x, p[[1]]), dim(x)))

matplot(res.2$time, res.2$mass.leaf[[1]], type="l", col="black", lty=1,
        xlab="Time", ylab="Leaf mass")

## Boo - these don't match *at all*.
matplot(res.2$time, res.2$mass.leaf[[1]], type="l", col="black", lty=1,
        xlab="Time", ylab="Leaf mass")
matlines(output$patch_age$age, output$strategies[[1]]$bound_m,
         type="l", lty=1,
         xlab="Time", ylab="Leaf mass", col="grey")
