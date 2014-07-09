library(tree)

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.0648406)))
p$add_strategy(new(Strategy, list(lma=0.1977910)))
p$seed_rain <- c(1.1, 1.1)               # Starting rain.
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster
p$set_control_parameters(equilibrium_verbose())

## Work out what the equilibrium seed rain is.  This takes a while!  If
## using the second set of seed_rain values this should converge
## quickly once the cohort schedule is worked out.
res <- equilibrium_seed_rain(p)

## Set the seed rain back into the parameters
p$seed_rain <- res$seed_rain[,"out"]

schedule <- res$schedule

## Resident lma values:
lma_res <- sapply(seq_len(p$size), function(i) p[[i]]$parameters$lma)

## OK, good to reasonable accuracy; the seed rain values should be
## close enough to 1.
cmp <- landscape("lma", lma_res, p, schedule)
cmp - 1 # This should be about zero, plus or minus 1e-5 or so.

# Mutant LMA values, in increasing numbers, to test how the time
# requirements scale with the number of strategies.  Should be
# sublinear.
lma_2  <- seq_log(0.03, 0.8, 2)
lma_4  <- seq_log(0.03, 0.8, 4)
lma_8  <- seq_log(0.03, 0.8, 8)
lma_16 <- seq_log(0.03, 0.8, 16)
lma_32 <- seq_log(0.03, 0.8, 32)
lma_64 <- seq_log(0.03, 0.8, 64)

(t_2  <- system.time(w_2  <- landscape("lma", lma_2,  p, schedule)))
(t_4  <- system.time(w_4  <- landscape("lma", lma_4,  p, schedule)))
(t_8  <- system.time(w_8  <- landscape("lma", lma_8,  p, schedule)))
(t_16 <- system.time(w_16 <- landscape("lma", lma_16, p, schedule)))
(t_32 <- system.time(w_32 <- landscape("lma", lma_32, p, schedule)))
(t_64 <- system.time(w_64 <- landscape("lma", lma_64, p, schedule)))

tt <- c(t_2[["elapsed"]],
        t_4[["elapsed"]],
        t_8[["elapsed"]],
        t_16[["elapsed"]],
        t_32[["elapsed"]],
        t_64[["elapsed"]])
n <- 2^seq_along(tt)

## This is largely OK, but the increased per-mutant time is
## concerning:
##+ time_per_mutant
plot(n, tt / n, log="x")

##+ fitness_landscape
plot(w_64 ~ lma_64, type="l", log="x")
abline(v=lma_res, lty=3, col="grey")

# Fitness of the mutants in an empty landscape:
w_empty <- landscape_empty(lma_16, p, schedule)

##+ fitness_landscape_empty
plot(w_empty ~ lma_16, log="x")
