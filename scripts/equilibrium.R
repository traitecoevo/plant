library(tree)
library(parallel)

run <- function(seed_rain_in, p, schedule) {
  p$seed_rain <- seed_rain_in
  run_ebt(p, schedule)$seed_rains
}

run_new_schedule <- function(w, p, schedule=NULL) {
  p$seed_rain <- w
  res <- build_schedule(p, schedule)
  unname(attr(res, "seed_rain")[,"out"])
}

## Try to establish what the equilubrium seed rain is.
p <- ebt_base_parameters()
p$add_strategy(new(Strategy, list(lma=0.08)))
p$seed_rain <- 1.1

## # 1: Approach to equilibrium:
res <- equilibrium_seed_rain(p)
seed_rain_eq <- unname(res[["seed_rain"]][,"in"])
schedule1 <- res[["schedule"]]$copy()

## Sanity and time check
system.time(delta <- run(seed_rain_eq, p, schedule1) -
            unname(res[["seed_rain"]][,"out"]))
delta

## Plot the approach:
approach <- t(sapply(attr(res, "progress"), "[[", "seed_rain"))

## From a distance, these both hone in nicely on the equilibrium, and
## rapidly, too.
##+ equilibrium_approach
r <- range(approach)
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach, pch=19, cex=.5, type="o")

## Zoom in on the last few points:
##+ equilibrium_approach_detail
r <- seed_rain_eq + c(-1, 1) * 0.03
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach, pch=19, cex=.5, type="o")

## Zoom in on the last few points to see where the insability kicks
## in:
##+ equilibrium_approach_wow_such_detail
r <- seed_rain_eq + c(-1, 1) * 0.000003
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach, pch=19, cex=.5, type="o")

## # 2: Near the equilibrium point:

## Then, in the vinicity of the root we should look at what the curve
## actually looks like, without adaptive refinement.
dr <- 2 # range of input to vary (plus and minus this many seeds)
seed_rain_in <- seq(seed_rain_eq - dr,
                    seed_rain_eq + dr, length=31)
seed_rain_out <- unlist(mclapply(seed_rain_in, run, p, schedule1))

fit <- lm(seed_rain_out ~ seed_rain_in)

## Here is input seeds vs output seeds:
##+ seeds_in_seeds_out
plot(seed_rain_in, seed_rain_out, xlab="Incoming seed rain",
     ylab="Outgoing seed rain", las=1)
abline(0, 1, lty=2, col="grey")
abline(fit, lty=2)
cobweb(approach)

## See instability.R for more exploration of this:
##+ seeds_in_seeds_out_instability
plot(seed_rain_in, resid(fit),
     xlab="Incoming seed rain", ylab="Residual seed rain",
     las=1, pch=1)
abline(h=0, v=seed_rain_eq)

## # 4. Global function shape
seed_rain_in_global <- seq(1, max(approach), length.out=51)

## This takes quite a while to compute.
##+ cache=TRUE
seed_rain_out_global <-
  unlist(mclapply(seed_rain_in_global, run_new_schedule, p))

## This is pretty patchy, which is due to incompletely refining the
## cohort schedule, I believe.  Tighten `schedule_eps` to make the
## curve smoother, at the cost of potentially a lot more effort.
##+ seeds_in_seeds_out_global
plot(seed_rain_in_global, seed_rain_out_global,
     las=1, type="l",
     xlab="Incoming seed rain", ylab="Outgoing seed rain")
abline(0, 1, lty=2, col="grey")
abline(fit,  lty=2)
cobweb(approach, lty=3)

## # 5. Multiple species at once:
p <- ebt_base_parameters()
p$add_strategy(new(Strategy, list(lma=0.08)))
p$add_strategy(new(Strategy, list(lma=0.15)))
p$seed_rain <- c(1.1, 1.1)               # Starting rain.

res <- equilibrium_seed_rain(p)
seed_rain_eq <- res[["seed_rain"]][,"out"]

progress_seed_rain <- lapply(attr(res, "progress"), "[[", "seed_rain")
approach <- lapply(seq_len(p$size), function(i)
                   t(sapply(progress_seed_rain, function(x) x[i,])))

## From a distance, these both hone in nicely on the equilibrium, and
## rapidly, too.
##+ approach_two_species
r <- range(unlist(approach))
plot(approach[[1]], type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cols <- c("black", "red")
for (i in seq_along(approach)) {
  cobweb(approach[[i]], pch=19, cex=.5, type="o", col=cols[[i]])
}
abline(v=seed_rain_eq, col=1:2, lty=3)
