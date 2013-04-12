source("helper-tree.R")
options(error=traceback)

context("Patch [CohortDiscrete]")

p <- new(Parameters)
p$add_strategy(new(Strategy))

## A plant that will be the same in terms of strategy (and initial
## mass).
cmp <- new(Plant, p$get_strategy(0))

patch.p <- new(Patch,  p)
patch.c <- new(PatchC, p)

expect_that(patch.p$height_max, is_identical_to(0.0))
expect_that(patch.c$height_max, is_identical_to(0.0))

## Add a single seed.
patch.p$add_seeds(1)
patch.c$add_seeds(1)

plants.p <- patch.p$get_plants()
plants.c <- patch.c$get_plants()

expect_that(length(plants.c), equals(1))
expect_that(length(plants.c[[1]]), equals(1))
expect_that(length(plants.c[[c(1,1)]]), equals(1))
expect_that(class(plants.c[[c(1,1)]]),
            equals("Rcpp_CohortDiscrete", check.attr=FALSE))
expect_that(class(plants.p[[c(1,1)]]),
            equals("Rcpp_Plant", check.attr=FALSE))

## Then clear both 
patch.p$clear()
patch.c$clear()

## Add a *pair* of seeds.
patch.p$add_seeds(2)
patch.c$add_seeds(2)

## And expect that the size of the system is 6 with no cohorts, and 3
## with cohorts.
expect_that(patch.p$ode_size, equals(6))
expect_that(patch.c$ode_size, equals(3))
## but the number of individuals is 2 in both cases:
expect_that(patch.p$n_individuals, equals(2))
expect_that(patch.c$n_individuals, equals(2))

## Now, grow this pair of seeds deterministically for a bit.
f <- function(patch, t) {
  res <- list()
  while ( patch$age < t ) {
    patch$step_deterministic()
    res <- c(res, list(c(patch$age, patch$ode_values)))
  }
  do.call(rbind, res)
}

res.p <- f(patch.p, 10)
res.c <- f(patch.c, 10)

expect_that(res.p[,1:4], equals(res.c[,1:4]))
expect_that(res.p[,5:7], equals(res.p[,2:4]))

## Add another 3 seeds.
patch.p$add_seeds(3)
patch.c$add_seeds(3)

expect_that(patch.p$n_individuals, equals(5))
expect_that(patch.c$n_individuals, equals(5))

## Clean up and start again.
patch.p$clear()
patch.c$clear()
patch.p$add_seeds(1)
patch.c$add_seeds(1)

## Function to run a patch up to a population size or age, and
## accumulate results.
run <- function(patch, n, age, res=NULL) {
  collect <- function()
    list(list(patch$age, patch$n_individuals, patch$get_mass_leaf(0)))
  if ( is.null(res) )
    res <- collect()
  while ( patch$n_individuals <= n && patch$age <= age ) {
    patch$step()
    res <- c(res, collect())
  }
  res
}

set.seed(1)
res.p <- run(patch.p, 1, 15)
set.seed(1)
res.c <- run(patch.c, 1, 15)

## Getting odd error here that goes away eventually.  Possible that we
## have a memory issue.
gc() # work around a bug with finalisers (not our fault?)
expect_that(res.p, equals(res.c))

## Continue on until the populations are quite large.  These will
## stochastically diverge once at least one cohort has more than one
## individual in it (see design.md).
set.seed(1)
res.p <- run(patch.p, 1000, 100, res.p)
set.seed(1)
res.c <- run(patch.c, 1000, 100, res.c)

nth <- function(x, n, strict=TRUE) if (strict) x[[n]] else x[n]
first  <- function(x, strict=TRUE) nth(x, 1, strict)
second <- function(x, strict=TRUE) nth(x, 2, strict)
third  <- function(x, strict=TRUE) nth(x, 3, strict)

t.p <- sapply(res.p, first)
t.c <- sapply(res.c, first)
n.p <- sapply(res.p, second)
n.c <- sapply(res.c, second)
m.p <- sapply(res.p, third)
m.c <- sapply(res.c, third)

m1.p <- sapply(m.p, first)
m2.p <- sapply(m.p, second, FALSE)
m1.c <- sapply(m.c, first)
m2.c <- sapply(m.c, second, FALSE)

## Check that results *are* the same up to the first cohort with >1
## individual.
i <- which(diff(n.c) > 1)[1]
expect_that(res.c[1:i], is_identical_to(res.c[1:i]))

## The first divergence may not happen for a while after this though.
j.p <- i:length(n.p)
j.c <- i:length(n.c)

## Check deterministic bits after the divergence (under assumption of
## no death in either of the two focal individuals -- currently true
## with this random seed and parameters, but possibly dangerous to
## rely on).
expect_that(m1.c[j.c],
            equals(spline(t.p[j.p], m1.p[j.p], xout=t.c[j.c])$y))
expect_that(m2.c[j.c],
            equals(spline(t.p[j.p], m2.p[j.p], xout=t.c[j.c])$y,
                   tolerance=1e-5)) # smaller plant, lower abs accuracy

## Then in a much more complicated fashion, check that the slope of
## the interpolated population size is not different than the line y =
## x (i.e., intercept 0, slope 1).
nn.p <- spline(t.p[j.p], n.p[j.p], xout=t.c[j.c])$y
nn.c <- n.c[j.c]
fit <- coef(summary(lm(log(nn.p) ~ log(nn.c))))
expect_that(fit[1,"Estimate"] + c(-2, 2) * fit[1,"Std. Error"],
            contains(0))
expect_that(fit[2,"Estimate"] + c(-2, 2) * fit[2,"Std. Error"],
            contains(1))

## This plots the results, which are interesting enough now.

## plot(m1.p ~ t.p, log="y", type="o", cex=.5, pch=19)
## lines(m2.p ~ t.p, lty=2)

## lines(m1.c ~ t.c, type="o", cex=.5, pch=19, col="#FF000077")
## lines(m2.c ~ t.c, lty=2, col="#FF000077")

## plot(n.p ~ t.p, type="s", log="y")
## lines(n.c ~ t.c, type="s", col="#FF000077")

## dt <- diff(t)
## plot(t[-length(t)], dt, type="s", log="y")
