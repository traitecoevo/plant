## Towards evolutionary assembly.
library(tree)
library(assertthat)
library(lhs)
library(mvtnorm)

source("assembly-fun.R")

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast.control())
t.max <- p0$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)

seed_rain0    <- 1e-3 # rain for new arrivals
seed_rain.eps <- 1e-4 # below which, consider species extinct

bounds <- rbind(lma=10^c(-2.5, 1.5),
                hmat=c(0.1,100))
colnames(bounds) <- c("lower", "upper")

## VCV that explores fraction p of the bounds.
p <- 0.001
vcv <- diag(2) * p * as.numeric(diff(t(log(bounds))))

mean.n.mutants    <- 2
mean.n.immigrants <- 1

new.phenotypes <- make.new.phenotypes(mean.n.mutants, vcv,
                                      mean.n.immigrants, bounds)
births <- make.births(new.phenotypes, seed_rain0, times0)
deaths <- make.deaths(seed_rain.eps)
run <- make.run(p0, t.max)

sys <- list(traits=cbind(lma=0.1978791, hmat=16.5958691),
            seed_rain=505.55407,
            times=list(times0))

## Here goes.  This will run indefinitely, so kill it once the
## community has stabilised.
res <- list(list(sys))
repeat {
  message(sprintf("*** Step %d:", length(res)))
  print(cbind(sys$traits, rain=sys$seed_rain))
  sys <- run(sys)
  sys <- deaths(sys)
  sys <- births(sys)
  res <- c(res, list(sys))
}
