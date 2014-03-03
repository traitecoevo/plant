library(tree)
source("landscape2-fun.R")

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.0648406, hmat=27)))
p$add_strategy(new(Strategy, list(lma=0.1977910, hmat=27)))
p$seed_rain <- c(1.1, 1.1)               # Starting rain.
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)
schedule0 <- schedule.from.times(times0, 2L)

# Work out what the equilibrium seed rain is.
res <- equilibrium.seed.rain(p, schedule0, 10,
                             build.args=list(verbose=TRUE),
                             progress=TRUE, verbose=TRUE)

# Before exploring the fitness landscape, put the ODE solver into
# non-adaptive mode by enforcing the same times as before.  That does
# require rerunning things...
schedule <- attr(res, "schedule")
ebt <- run.ebt(p, schedule$copy())
schedule$ode_times <- ebt$ode_times

# Resident lma values:
lma.res <- sapply(seq_len(p$size), function(i) p[[i]]$parameters$lma)

# Some mutants:
lma <- seq.log(0.03, 0.8, 12)

# Fitness of the mutants (in the presence of the residents)
w <- landscape(lma, p, schedule)

plot(w ~ lma, log="x")
abline(v=lma.res, lty=3, col="grey")

# Fitness of the mutants in an empty landscape:
w.empty <- landscape.empty(lma, p, schedule)

plot(w.empty ~ lma, log="x")
