library(tree)
source("landscape2-fun.R")

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.0648406)))
p$add_strategy(new(Strategy, list(lma=0.1977910)))
p$seed_rain <- c(1.1, 1.1)               # Starting rain.
p$seed_rain <- c(348.738300488797, 3.07449195531246) # or start close
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)
schedule0 <- schedule.from.times(times0, 2L)

# Work out what the equilibrium seed rain is.  This takes a while!  If
# using the second set of seed_rain values this should converge
# quickly once the cohort schedule is worked out.
res <- equilibrium.seed.rain(p, schedule0, 10,
                             build.args=list(verbose=TRUE),
                             progress=TRUE, verbose=TRUE)
# Set the seed rain back into the parameters
p$seed_rain <- res$seed_rain[,"out"]

# And rerun again to get the ODE times.
schedule <- res$schedule
ebt <- run.ebt(p, schedule$copy())
schedule$ode_times <- ebt$ode_times

# Resident lma values:
lma.res <- sapply(seq_len(p$size), function(i) p[[i]]$parameters$lma)

## OK, good to reasonable accuracy; the seed rain values should be
## close enough to 1.
cmp <- landscape(lma.res, p, schedule)
cmp - 1 # This should be about zero, plus or minus 1e-5 or so.

# Mutant LMA values, in increasing numbers, to test how the time
# requirements scale with the number of strategies.  Should be
# sublinear.
lma.2  <- seq.log(0.03, 0.8, 2)
lma.4  <- seq.log(0.03, 0.8, 4)
lma.8  <- seq.log(0.03, 0.8, 8)
lma.16 <- seq.log(0.03, 0.8, 16)
lma.32 <- seq.log(0.03, 0.8, 32)
lma.64 <- seq.log(0.03, 0.8, 64)

(t.2  <- system.time(w.2  <- landscape(lma.2,  p, schedule)))
(t.4  <- system.time(w.4  <- landscape(lma.4,  p, schedule)))
(t.8  <- system.time(w.8  <- landscape(lma.8,  p, schedule)))
(t.16 <- system.time(w.16 <- landscape(lma.16, p, schedule)))
(t.32 <- system.time(w.32 <- landscape(lma.32, p, schedule)))
(t.64 <- system.time(w.64 <- landscape(lma.64, p, schedule)))

tt <- c(t.2[["elapsed"]],
        t.4[["elapsed"]],
        t.8[["elapsed"]],
        t.16[["elapsed"]],
        t.32[["elapsed"]],
        t.64[["elapsed"]])
n <- 2^seq_along(tt)
plot(n, tt / n, log="x")

plot(w.64 ~ lma.64, type="l", log="x")
abline(v=lma.res, lty=3, col="grey")

# Fitness of the mutants in an empty landscape:
w.empty <- landscape.empty(lma.16, p, schedule)

plot(w.empty ~ lma.16, log="x")
