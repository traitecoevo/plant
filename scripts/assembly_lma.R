library(tree)
library(assertthat)
library(lhs)
library(mvtnorm)

source("scripts/assembly-fun.R")

# Paremeters that vary among runs

assemble <- function(time.disturbance=10, prod=1.0, n.steps=1000, output.dir="output",
  p=0.001, mean.n.mutants=1, mean.n.immigrants=1){

filename <- sprintf("%s/assembly-%s-%s.rds", output.dir, time.disturbance, prod)
if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)

## Default strategy with no reproduction
strategy = new(Strategy)
strategy$set_parameters(structure(list(0.5, 0, strategy$parameters[["c_p1"]]*prod), names=c("c_r1", "c_r2", "c_p1")))

p0 <- new(Parameters)
p0$set_parameters(list(patch_area=1.0))
p0$set_control_parameters(fast.control())
p0$disturbance <- new(Disturbance, time.disturbance)

t.max <- p0$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)

seed_rain0    <- 1e-3 # rain for new arrivals
seed_rain.eps <- 1e-4 # below which, consider species extinct

bounds <- rbind(lma=10^c(-2.5, 1.5))

colnames(bounds) <- c("lower", "upper")

## VCV that explores fraction p of the bounds.
vcv <- diag(nrow(bounds)) * p * as.numeric(diff(t(log(bounds))))

new.phenotypes <- make.new.phenotypes(mean.n.mutants, vcv,
                                      mean.n.immigrants, bounds)
births <- make.births(new.phenotypes, seed_rain0, times0)
deaths <- make.deaths(seed_rain.eps)
run <- make.run(p0, t.max, strategy=strategy)

sys <- list(traits=cbind(lma=0.1978791),
            seed_rain=505.55407,
            times=list(times0))

## Take a given number of steps, saving output every so often
res <- list(list(sys))
for(i in 1:n.steps){
  message(sprintf("*** Step %d:", length(res)))
  print(cbind(sys$traits, rain=sys$seed_rain))
  sys <- run(sys)
  sys <- deaths(sys)
  sys <- births(sys)
  res <- c(res, list(sys))
  if(i %% 10)
	 saveRDS(res, filename)
}

saveRDS(res, filename)
}
