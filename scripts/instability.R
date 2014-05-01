library(tree)
library(parallel)

## Try to establish what the equilubrium seed rain is.
p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- 1.1                       # Starting rain.
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)
schedule0 <- schedule.from.times(times0)

run <- function(seed_rain.in, p, schedule) {
  p$seed_rain <- seed_rain.in
  run.ebt(p, schedule)$fitnesses
}

## TODO: This could save the ode times at the same time here.
run.collect <- function(seed_rain.in, p, schedule) {
  p$seed_rain <- seed_rain.in
  run.ebt.collect(p, schedule)
}

run.new.schedule <- function(w, p, schedule, build.args=list()) {
  p$seed_rain <- w
  build.args <- modifyList(list(nsteps=20, eps=1e-3, verbose=FALSE),
                           build.args)
  res <- build.schedule(p, schedule, build.args$nsteps, build.args$eps,
                        progress=FALSE, verbose=build.args$verbose)
  unname(attr(res, "seed_rain")[,"out"])
}

## # 1: Approach to equilibrium:
res <- equilibrium.seed.rain(p, schedule0, 10, progress=TRUE,
                             build.args=list(verbose=TRUE))
w.hat <- unname(res[["seed_rain"]][,"in"])
w.hat.out <- unname(res[["seed_rain"]][,"out"])
schedule1 <- res[["schedule"]]$copy()

approach <- t(sapply(attr(res, "progress"), "[[", "seed_rain"))

## Sanity and time check
system.time(w.cmp <- run(w.hat, p, schedule1)) # ~10s
w.cmp - w.hat.out # should be exactly zero

## Check that running the ebt with collection does not change things:
system.time(res.cmp <- run.collect(w.hat, p, schedule1)) # ~13s
res.cmp$fitness - w.hat.out

## Then, in the vinicity of the root we should look at what the curve
## actually looks like, without adaptive refinement.
dw <- 1 # range of input to vary (plus and minus this many seeds)
seed_rain.in  <- seq(w.hat - dw, w.hat + dw, length=31)
seed_rain.out <- unlist(mclapply(seed_rain.in, run, p, schedule1))
fit <- lm(seed_rain.out ~ seed_rain.in)

p2 <- p$copy() # modify fast control
p2$set_control_parameters(list(environment_light_tol=1e-6))
system.time(w.cmp2 <- run(w.hat, p2, schedule1)) # ~14s
seed_rain.out2 <- unlist(mclapply(seed_rain.in, run, p2, schedule1))
fit2 <- lm(seed_rain.out2 ~ seed_rain.in)

p3 <- p$copy() # modify fast control
p3$set_control_parameters(list(ode_tol_rel=1e-5, ode_tol_abs=1e-5))
system.time(w.cmp3 <- run(w.hat, p3, schedule1)) # ~26s
seed_rain.out3 <- unlist(mclapply(seed_rain.in, run, p3, schedule1))
fit3 <- lm(seed_rain.out3 ~ seed_rain.in)

p4 <- p$copy() # modify fast control
p4$set_control_parameters(list(plant_assimilation_rule=61))
system.time(w.cmp4 <- run(w.hat, p4, schedule1)) # ~17s
seed_rain.out4 <- unlist(mclapply(seed_rain.in, run, p4, schedule1))
fit4 <- lm(seed_rain.out4 ~ seed_rain.in)

p5 <- p2$copy() # modify high-precision light environment
p5$set_control_parameters(list(plant_assimilation_rule=61))
system.time(w.cmp5 <- run(w.hat, p5, schedule1)) # ~20s
seed_rain.out5 <- unlist(mclapply(seed_rain.in, run, p5, schedule1))
fit5 <- lm(seed_rain.out5 ~ seed_rain.in)

# fast, better light environment, better ode, better assimilation,
# better assimilation and better light environment
cols <- c("black", "red", "blue", "green4", "purple")

## Here is input seeds vs output seeds:
##+ seeds_in_seeds_out
matplot(seed_rain.in,
        cbind(seed_rain.out,  seed_rain.out2, seed_rain.out3,
              seed_rain.out4, seed_rain.out5),
        col=cols, las=1, pch=1,
        xlab="Incoming seed rain", ylab="Outgoing seed rain")
abline(0, 1, lty=2, col="grey")
abline(fit, lty=2, col=cols[1])
abline(fit2, lty=2, col=cols[2])
abline(fit3, lty=2, col=cols[3])
abline(fit4, lty=2, col=cols[4])
abline(fit5, lty=2, col=cols[5])

## The instability is easier to see here.  There are two patterns; the
## blue and black lines (using the low order integration) have a
## switch point.  The other cases have less systematic error, though
## general error on the same order.  Oddly increasing both the
## assimilation integration order and the light environment
## sensitivity did worse than just increasing the light environment
## sensitivity.
##+ seeds_in_seeds_out_residual
matplot(seed_rain.in,
        cbind(resid(fit), resid(fit2), resid(fit3), resid(fit4), resid(fit5)),
        xlab="Incoming seed rain", ylab="Residual seed rain",
        las=1, pch=1, col=cols, type="o", lty=1, cex=.5)
abline(h=0, v=w.hat)

## I wonder if the issue is simply the trapezium integration at the
## end?  Plot the components of fitness for the two points that are
## close together and see where that difference is coming from --
## through out, just a single point, or something more subtle.  We can
## either get this from the code in test-ebt.R or using
## EBT::fitness_cohort (for which we'll need to tweak ebt.collect).

## It's easy to grab the points that differ most:
##+ seeds_in_seeds_out_biggest_difference
i <- which.max(abs(diff(seed_rain.out)))
j <- i + 1
matplot(seed_rain.in,
        cbind(resid(fit), resid(fit2), resid(fit3), resid(fit4), resid(fit5)),
        xlab="Incoming seed rain", ylab="Residual seed rain",
        las=1, pch=1, col=cols, type="o", lty=1, cex=.5)
abline(h=0, v=w.hat)
points(seed_rain.in[c(i,j)], resid(fit)[c(i, j)], pch=19)

res.i <- run.collect(seed_rain.in[[i]], p, schedule1)
res.j <- run.collect(seed_rain.in[[j]], p, schedule1)

## Quick check -- should be zero:
seed_rain.out[c(i, j)] - c(res.i$fitness, res.j$fitness)

## These have to the be the same:
identical(res.i$time, res.j$time)
identical(dim(res.i$species[[1]]), dim(res.j$species[[1]]))

## Just one difference in the light environment.
n.light <- cbind(sapply(res.i$light.env, nrow),
                 sapply(res.j$light.env, nrow))
sum(n.light[,1] != n.light[,2])

## But is that important?  How many in neighbouring ones too?
res.i1 <- run.collect(seed_rain.in[[i-1]], p, schedule1)
res.j1 <- run.collect(seed_rain.in[[j+1]], p, schedule1)
n.light1 <- cbind(sapply(res.i1$light.env, nrow),
                  n.light,
                  sapply(res.j1$light.env, nrow))
## OK, looks like it could be important:
rowSums(diff(t(n.light1)))
