## Adaptive integration causes problems with the smoothness of the
## seed rain in/seed rain out function.
library(tree)
library(parallel)

## Try to establish what the equilibrium seed rain is.  First, a
## non-adaptive set of parameters
p_n <- new(Parameters)
p_n$add_strategy(new(Strategy))
p_n$seed_rain <- 1.1                       # Starting rain.
p_n$set_parameters(list(patch_area=1.0))   # See issue #13
p_n$set_control_parameters(fast.control()) # A bit faster
p_n$set_control_parameters(equilibrium_verbose())

## And a set that is identical except for the adaptive assimilation is
## turned on:
p_a <- p_n$copy()
p_a$set_control_parameters(list(plant_assimilation_adaptive=TRUE))

run <- function(seed_rain_in, p, schedule) {
  p$seed_rain <- seed_rain_in
  run.ebt(p, schedule)$fitnesses
}

run_new_schedule <- function(w, p, schedule=NULL) {
  p$seed_rain <- w
  res <- build_schedule(p, schedule)
  unname(attr(res, "seed_rain")[,"out"])
}

## # 1: Approach to equilibrium:
res_n <- equilibrium_seed_rain(p_n)
res_a <- equilibrium_seed_rain(p_a)

approach_n <- t(sapply(attr(res_n, "progress"), "[[", "seed_rain"))
approach_a <- t(sapply(attr(res_a, "progress"), "[[", "seed_rain"))

cols <- c(n="black", a="red")
w_hat <- unname(res_n[["seed_rain"]][,"out"])

## From a distance, these both hone in nicely on the equilibrium, and
## rapidly, too.
##+ approach
r <- range(approach_n, approach_a)
plot(approach_n, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach_n, col=cols[["n"]], pch=19, cex=.5, type="o")
cobweb(approach_a, col=cols[["a"]], pch=19, cex=.5, type="o")

## Here is the painful lack of convergence for the
## with-adaptive-integration case (red); it's converging to a *region*
## but moving essentially stochastically within that region.
##+ approach_detail
r <- w_hat + c(-1, 1) * 0.5
plot(approach_n, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach_n, col=cols[["n"]], pch=19, cex=.5, type="o")
cobweb(approach_a, col=cols[["a"]], pch=19, cex=.5, type="o")

## # 2: Near the equilibrium point:

## Then, in the vinicity of the root we should look at what the curve
## actually looks like, without adaptive refinement.
dw <- 2 # range of input to vary (plus and minus this many seeds)
seed_rain_in  <- seq(w_hat - dw, w_hat + dw, length=31)

schedule1 <- res_n$schedule
seed_rain_out_n <- unlist(mclapply(seed_rain_in, run, p_n, schedule1))
seed_rain_out_a <- unlist(mclapply(seed_rain_in, run, p_a, schedule1))

fit_n <- lm(seed_rain_out_n ~ seed_rain_in)
fit_a <- lm(seed_rain_out_a ~ seed_rain_in)

##+ seeds_in_seeds_out
matplot(seed_rain_in, cbind(seed_rain_out_n, seed_rain_out_a),
        xlab="Incoming seed rain", ylab="Outgoing seed rain",
        las=1, pch=1, col=cols)
abline(0, 1, lty=2, col="grey")
abline(fit_n, lty=2, col=cols[["n"]])
abline(fit_a, lty=2, col=cols[["a"]])
cobweb(approach_n, col=cols[["n"]])
cobweb(approach_a, col=cols[["a"]])

## # 3: Global function shape
seed_rain_in_global <- seq(1, w_hat + 10, length=51)

seed_rain_out_global_n <-
  unlist(mclapply(seed_rain_in_global, run_new_schedule, p_n))
seed_rain_out_global_a <-
  unlist(mclapply(seed_rain_in_global, run_new_schedule, p_a))

##+ approach_global
matplot(seed_rain_in_global,
        cbind(seed_rain_out_global_n, seed_rain_out_global_a),
        las=1, type="l", col=cols, lty=1,
        xlab="Incoming seed rain", ylab="Outgoing seed rain")
abline(0, 1, lty=2, col="grey")
abline(fit_n, lty=2, col=cols[["n"]])
abline(fit_a, lty=2, col=cols[["a"]])
cobweb(approach_n, col=cols[["n"]])
cobweb(approach_a, col=cols[["a"]])
