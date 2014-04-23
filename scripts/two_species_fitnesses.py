import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
import numpy as np

# Load tree in R
robjects.r('source("~/.Rprofile")') # hopefully this sets .libPaths
robjects.r("library(tree)")

def run_simmo(pop_lma_arr, seed_rain_arr = None):
    """ Run tree to calculate a fitness landscape with given species.

        Based on landscape2.R.
    """
    robjects.r("p <- new(Parameters)")

    # Starting seed rain
    if seed_rain_arr == None:
        seed_rain_arr = 1.1 * np.ones_like(pop_lma_arr)

    assert len(pop_lma_arr) == len(seed_rain_arr)

    for lma in pop_lma_arr:
        # Add species into simulation
        robjects.r("p$add_strategy(new(Strategy, list(lma=%f)))" % lma)

    robjects.r("p$seed_rain <- c(%s)" % str(list(seed_rain_arr))[1:-1])

    robjects.r("""
        p$set_parameters(list(patch_area=1.0))   # See issue #13
        p$set_control_parameters(fast.control()) # A bit faster

        t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
        times0 <- cohort.introduction.times(t.max)
        """)
    robjects.r("schedule0 <- schedule.from.times(times0, %iL)" % len(pop_lma_arr))
    robjects.r("""
        res <- equilibrium.seed.rain(p, schedule0, 10,
                                     build.args=list(verbose=TRUE),
                                     progress=TRUE, verbose=TRUE)
        # Set the seed rain back into the parameters
        p$seed_rain <- res$seed_rain[,"out"]

        # And rerun again to get the ODE times.
        schedule <- res$schedule
        ebt <- run.ebt(p, schedule$copy())
        schedule$ode_times <- ebt$ode_times

        # Mutant LMA values, in increasing numbers, to test how the time
        # requirements scale with the number of strategies.  Should be
        # sublinear.
        lma.f <- seq_log(0.03, 0.8, 64)

        w.f <- landscape(lma.f, p, schedule)
        """)

    out_seed_rain = np.array(robjects.r('res$seed_rain[,"out"]'))
    fitness = np.array(robjects.r['w.f'])
    lma = np.array(robjects.r['lma.f'])
    return out_seed_rain, fitness, lma


seed_rains = {}
fitnesses = {}

last_seed_rain = [1.1, 1.1]
fixed_lma = 0.08

for lma in np.logspace(np.log(0.03), np.log(0.8), 64, base=np.e)[5:-5]:
    # Avoid seed rains getting stuck on zero
    input_seed_rain = [max(1.1, s) for s in last_seed_rain]

    last_seed_rain, fitness, _ = run_simmo([lma, fixed_lma],
                                            input_seed_rain)
    seed_rains[lma] = last_seed_rain
    fitnesses[lma] = fitness


## Assorted sloppy bits that I type into IPython --pylab to run this:

# %run /path/to/this/script.py

## Save raw results to avoid recomputation (nb hard-coded filenames)
# import pickle
# with open("~/Data/output/tree/two_species_008_lmas.pickle", 'w') as f:
#     pickle.dump(np.array(robjects.r['lma.f']), f)
# with open("~/Data/output/tree/two_species_008_fitnesses.pickle", 'w') as f:
#     pickle.dump(fitnesses, f)
# with open("~/Data/output/tree/two_species_008_seed_rains.pickle", 'w') as f:
#     pickle.dump(seed_rains, f)

# # Prepare results and plot
# all_fit = np.logspace(np.log(0.03), np.log(0.8), 64, base=np.e)
# fs = np.zeros((len(fitnesses), len(fitnesses.values()[0])))
# lmas = seed_rains.keys()
# for i, k in enumerate(sort(fitnesses.keys())):
#     fs[i] = fitnesses[k]
# figure()
# loglog()
# scatter(lmas, lmas, [20*log10(seed_rains[l][0]) for l in lmas], 'w', edgecolors='none')
# scatter(fixed_lma*np.ones_like(lmas), lmas, [20*log10(seed_rains[l][1]) for l in lmas], 'w', edgecolors='none')
# imshow((fs-1)/np.max(fs, axis=1, keepdims=True), extent=(0.03, 0.8, max(fitnesses.keys()), min(fitnesses.keys())))
# colorbar()
