import tree

time_disturbance = 11.3
slope = 1.3
lma = [0.018, 0.09485587]
seed_rain = [450.44073695145, 1.53529930816904]

m = tree.TreeModel(time_disturbance, slope)

# First, compute fitness for a given set of species traits and
# densities, but *not* at equilibrium:
m.set_residents(lma, seed_rain, False)

# This computes the true fitness values at the vector lma:
w = m.fitness(lma)

# This will construct an approximate landscape.  By default this will
# use our GP approach, which samples a bunch of points sequentially.
# It takes a while and is fairly chatty.  We can turn it down if that
# gets annoying.
w = m.fitness_approximate(lma)

# Generate a series of 50 points to estimate fitness for:
lma2 = m.grid(50)
# And compute the fitness
w2 = m.fitness_approximate(lma2)

# At this point I get baffled because in R I'd want to plot these, but
# plotting in Python is a black art to me :)

# An object representing the GP is also available:
gp = m.gp
gp.X # lma locations (but logged)
gp.y # fitness values at these locations
gp.predict(0) # predict values for a given log(lma) - returns (mean, variance)

# Bear in mind that is a python interface to an R interface to some
# C++ code so may not be the flashest way of working with this.
# But if you need specific things from the GP object, then I can
# expose them in both places and you can get access easily enough.

# For completeness, this is the version where we stabilise the species
# abundances to get demographic equilibrium.  This is what we mostly
# use, but it sounded like you were more interested in the version
# where the abundances were free.
m.set_residents(lma, seed_rain, True)
# Force generation of the approximate landscape:
m.setup_fitness_approximate()
# As before:
gp = m.gp
gp.predict(0)
m.fitness_approximate([1])
