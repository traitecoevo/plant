library(tree)
library(parallel)

## With the current model/set of parameters, things here are nowhere
## near as bad as they used to be.  However, there is still some
## useful information here.

## Try to establish what the equilubrium seed rain is.
p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- 1.1                       # Starting rain.
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster
p$set_control_parameters(equilibrium_verbose())

run <- function(seed_rain_in, p, schedule) {
  p$seed_rain <- seed_rain_in
  run.ebt(p, schedule)$fitnesses
}

run_collect <- function(seed_rain_in, p, schedule) {
  p$seed_rain <- seed_rain_in
  run.ebt.collect(p, schedule)
}

run_new_schedule <- function(w, p, schedule=NULL) {
  p$seed_rain <- w
  res <- build_schedule(p, schedule)
  unname(attr(res, "seed_rain")[,"out"])
}

## # 1: Approach to equilibrium:
res <- equilibrium_seed_rain(p)
w_hat <- unname(res[["seed_rain"]][,"in"])
w_hat_out <- unname(res[["seed_rain"]][,"out"])
schedule1 <- res[["schedule"]]$copy()

approach <- t(sapply(attr(res, "progress"), "[[", "seed_rain"))

## Sanity and time check
system.time(w_cmp <- run(w_hat, p, schedule1)) # ~10s
w_cmp - w_hat_out # should be exactly zero

## Check that running the ebt with collection does not change things:
system.time(res_cmp <- run_collect(w_hat, p, schedule1)) # ~13s
res_cmp$fitness - w_hat_out

## Then, in the vinicity of the root we should look at what the curve
## actually looks like, without adaptive refinement.
dw <- 1 # range of input to vary (plus and minus this many seeds)
seed_rain_in  <- seq(w_hat - dw, w_hat + dw, length=31)
seed_rain_out <- unlist(mclapply(seed_rain_in, run, p, schedule1))
fit <- lm(seed_rain_out ~ seed_rain_in)

p2 <- p$copy() # modify fast control
p2$set_control_parameters(list(environment_light_tol=1e-6))
system.time(w_cmp2 <- run(w_hat, p2, schedule1)) # ~20s
seed_rain_out2 <- unlist(mclapply(seed_rain_in, run, p2, schedule1))
fit2 <- lm(seed_rain_out2 ~ seed_rain_in)

p3 <- p$copy() # modify fast control
p3$set_control_parameters(list(ode_tol_rel=1e-5, ode_tol_abs=1e-5))
system.time(w_cmp3 <- run(w_hat, p3, schedule1)) # ~21s
seed_rain_out3 <- unlist(mclapply(seed_rain_in, run, p3, schedule1))
fit3 <- lm(seed_rain_out3 ~ seed_rain_in)

p4 <- p$copy() # modify fast control
p4$set_control_parameters(list(plant_assimilation_rule=61))
system.time(w_cmp4 <- run(w_hat, p4, schedule1)) # ~13s
seed_rain_out4 <- unlist(mclapply(seed_rain_in, run, p4, schedule1))
fit4 <- lm(seed_rain_out4 ~ seed_rain_in)

p5 <- p2$copy() # modify high-precision light environment
p5$set_control_parameters(list(plant_assimilation_rule=61))
system.time(w_cmp5 <- run(w_hat, p5, schedule1)) # ~24s
seed_rain_out5 <- unlist(mclapply(seed_rain_in, run, p5, schedule1))
fit5 <- lm(seed_rain_out5 ~ seed_rain_in)

# fast, better light environment, better ode, better assimilation,
# better assimilation and better light environment
cols <- c("black", "red", "blue", "green4", "purple")

## Here is input seeds vs output seeds:
##+ seeds_in_seeds_out
matplot(seed_rain_in,
        cbind(seed_rain_out,  seed_rain_out2, seed_rain_out3,
              seed_rain_out4, seed_rain_out5),
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
matplot(seed_rain_in,
        cbind(resid(fit), resid(fit2), resid(fit3), resid(fit4), resid(fit5)),
        xlab="Incoming seed rain", ylab="Residual seed rain",
        las=1, pch=1, col=cols, type="o", lty=1, cex=.5)
abline(h=0, v=w_hat)
