# Overall

This represents our "version 1" of the model.  In it we will have a
deterministic model and a stochastic individual based simulation.
However, the scope is very much the same as the original "EBT" code.

 - no explit space at the within patch level (but possibly at the
   between patch level).
 - one size dimension (but possibly more for the individual based
   model).
 - large number of possibly varying traits (in contrast to the 4
   allowed in the original EBT implementation).

I've tried to abstract things to a useful degree, but not guess at how
much further abstraction would lead to more flexible modelling in the
future.  In particular, I don't think that the EBT code to control the
elevator-boxcar-train is sufficiently general to throw at other
problems as-is, and would certainly require further modification.

# General structure of components

For the stochastic model, we have:

* A `Metacommunity`, which can be advanced through time, containing...
* ...a number of `Patch`es.  Each of which has...
  - ...an `Environment`, containing information on patch survival,
    incoming seed rain, and the light environment.
  - ...a number of `Species`, each of which has...
* ...a number of `Plant`s.  The plant requires many parameters to
  define it's growth strategy, defined in a...
* ...`Strategy` object, which is unique per-`Species`.

For the EBT, the structure is very similar, except that at the top
level we have an `EBT` object that can be advanced through time that
contains a *single* `Patch` object, representing an idealised patch.

There are a number of different types of `Plant`; there is the plain
`Plant` that defines the basic approach, a `CohortDiscrete` that
contains an optimisation for the stochastic model, and there is a
`CohortTop` that works with the EBT.  In addition, there are a couple
of experimental plant types that don't yet work (`PlantApprox` and
`PlantSpline` -- just ignore these).  There will eventually be a
`CohortMean` that works with the EBT, but with a different
approximation.

# Description of model components, from bottom up

These are all in namespace `model`.

## Individuals (Plant, etc)

The "individuals" in the individual based model are very similar to
cohorts in the deterministic model; they have a state (mass, leaf
area, etc) and are capable of computing their physiological rates
given a light environment.

### Key variables

The key *variables* that define an individual are

* **Plant size**: Currently the only size variable is "height".  Note
  that this differs from the paper where "leaf mass" was the key size
  variable.  See `details.md` for information about consequences of
  this change of variables.
* **Mortality**: We don't directly track the probability of mortality,
  but instead track the rate of a Poisson process describing
  mortality.  So if $m$ is this rate variable, then the probability of
  dying is given as $1 - \exp(-m)$, and survival as $\exp(-m)$.
* **Fecundity**: We track the number of seeds produced, as a decimal
  number.  So a plant might have generated enough reproductive biomass
  for '12.4' seeds, rather than just integer numbers of seeds.
  
These variables are retrieved by the methods `height`, `mortality` and
`fecundity` and set by `set_height`, `set_mortality` and
`set_fecundity`, respectively.

### Impact of Plant on external environment

We need to be able to interrogate plants to build up a picture of
their impact on the light environment.  The only thing we need here is
the amount of leaf area deployed above a given height (via method
`leaf_area_above`).  The `Patch` will compute leaf area index given
this information.  Note that this would completely change in a spatial
model.

### Impact of external environment on Plant

There are two (related) places where the plant is affected by the
environment (principly the light environment).

First, in computing the rates of change of all key variables.  This is
done through the `compute_vars_phys` method, which takes an
`Environment` object.  rates of change of the three key variables
(`height`, `mortality`, `fecundity`) are computed.  The last set
values can be accessed by the methods `height_rate`, `mortality_rate`
and `fecundity_rate`, which gives the *rate of change with respect to
time* of each key variable.  This is a bit confusing for `fecundity`,
which is a rate itself (seeds/time), so what is returned is the
$d/dt[seeds/t]$.

### ODE interface

There is a common ODE stepping interface (`ode::OdeTarget`) that
`Plant` objects implement.  Plants implement the methods

* `ode_size`: return number of underlying variables (3)
* `ode_values_set`: set underlying variables
* `ode_values`: return underlying variables
* `ode_rates`: return rates of change of underlying variables

This is explained in the ODE section of this document.

### Individual based model

There are methods `births` and `deaths` for querying how many seeds
have been produced and checking if the individual has died.  These
reset the variables `fecundity` and `mortality` as they are used.

Because `fecundity` can indicate fractional numbers of seeds, `births`
returns the integral part, leaving the fractional part for future time
steps.

The method `died` looks to see if the plant has died since death was
last checked (so `mortality` increases over this period, then we run a
Bernoulli trial with probability `1-exp(-mortality)` to test for
death).  Then we reset `mortality` back to zero.

### EBT/Cohorts

*[TODO] handling of `density` and `seeds_survival_weighted` is still
really unclear and may change.  Going to need to push these into new
variables, I think.*

In addition to the three key variables, there is a fourth:

* `density`: related to the density (with respect to area) of
  individuals of this cohort size in the patch.  The *actual* density
  is given as `exp(-density)`
  
This is retreived with the `density` method, set with `set_density`
and `compute_vars_phys` sets the rate of change, which can be accessed
with `density_rate`.  Because this is an extra variable, `ode_size`
for `CohortTop` returns 4, rather than 3.

The `fecundity` variable is also swapped out for a variable
`seeds_survival_weighted`.  This is the lifetime production of seeds
of a plant (like `fecundity`) but weighted by survival.  There is some
documentation work needed here to explain how this works (TODO).

In contrast with the individual-based `Plant`, a `CohortTop`:

* `height` increases simply with the height growth rate (as in `Plant`).
* `mortality` increases with the mortality rate (as in `Plant`)
* `density` increases at rate $\partial growth_rate/ \partial m + death_rate$

Initial conditions do differ (relative to `Plant`):

* `height` unchanged (height of freshly germinated seed)
* `mortality` starts at `exp(1 - germination probability)`
* `density` starts at ?
* `seeds_survival_weighted` unchanged (0)

[TODO] I think that `seeds_survival_weighted` would be best to put
into the fecundity variable itself; it's not clear why it's not
currently there, actually...

The leaf area calculation is quite different in the cohort calculation
to the individual based simulation.  In fact, we're not actually sure
that it's correct yet (TODO).

### CohortDiscrete

This is a simple optimisation for the individual-based model.  If
multiple individuals are born at the same time, then we only have to
compute their properties once (growth is deterministic).  This simply
adds a constructor that takes an integer argument indicating how many
individals are in a cohort, and modifies a few methods where the
output depends on the number of individuals in the cohort.  The number
can also be set using the method `set_n_individuals`.

This trick affects the interpretation of a few methods (though this is
fairly transparent to classes that use `CohortTop` rather than
`Plant`) 

* `leaf_area_above`: returns `n_individuals * leaf_area_above`, as
  there are `n_individuals` at the given height.
* `offspring`: similarly, returns `n_individuals * offspring`
* `died`: only returns `true` (for death) when *all* individuals have
  died.  However, the cohort will decrease in size over time (changing
  `n_individuals`) as individuals die.
  
### Parameter handling via `Strategy`

The class `Strategy` defines the (many) parameters for the model.
This class has no public data members, but has `Plant` as a friend (so
`Plant` can access all of the members.  I chose that approach because
`Strategy` and `Plant` are very tightly linked (changing one can
easily require changing the other).  At the same time, nothing else
needs to know anything about the parameters that define a plant
strategy.

All parameters have default values following those in the paper.  The
four "core" traits have default values added too; these came from the
orginal EBT implementation (rationalle is: it doesn't really make much
sense to make a strategy that doesn't know its lma).

### Other notes

At the moment, all types are implemented by inheriting from Plant, as
we need most of the machinery there.  However, there are some aspects
of `Plant` that are tied up with the stochastic version (`offspring`
and `died` in particular) that aren't needed in the EBT version.
However, `Species` will also have these, so need to be careful about
how these are exposed.

[TODO] It could be that a base plant type that doesn't know any of
these things would be good to use here.  Because of the way that
`Species` and `Patch` are implemented (via templating), it's OK if
these methods are not included.  At present this is enforced by
runtime errors via `Rf_error`.

## Species

Above the level of species, there are two places that "species"
appears as a concept:

1. as a structure/class of parameters, corresponding to a different
   strategy.
2. as a book-keeping device, collecting individuals within a
   population

Different `Strategy` instances correspond to different conceptual
species.  The class `Species` is used for book-keeping and collects up
a number of individuals (`Plant`, `CohortDiscrete`, or `CohortTop`).

Depending on what a species is composed of, some methods should be
handled differently -- for example, `leaf_area_above`, `births`,
`deaths`.  Some methods may not even make sense.  At the moment these
are done via template specialisation.

The method `size()` returns the number of individuals contained (not
necessarily plants, as both `CohortDiscrete` and `CohortTop` stand in
for potentially more than one plant).

The main reason for having things in species are:

1. It makes it easy to see if different types are successful or not,
   because we only care about things like seed production at the level
   of species.
2. We can keep the species sorted by size and take advantage of that
   in a few places (e.g., working out the height of the largest seed,
   or looking only at tall plants when characterising the light
   environment).
3. The EBT approach really requires this level of grouping (because
   cohorts are representatives of a species).  This may weaken for the
   stochastic model where mutation is very common, as every mutation
   will create a new "species".
   
Most of the methods of `Species` either aggregate/compute quantities
over the set of plants (e.g., `height_max` returns the height of the
tallest individual, `leaf_area_above` computes the leaf area above a
hight over all individuals) or iterates over individuals (e.g.,
`compute_vars_phys`, which just runs the method with the same name
over all individuals within the species).  The methods `births` and
`deaths` loop over all plants checking for births and deaths,
returning the total number of births/deaths over the population.
(NOTE: the semantics here are odd for `Species<CohortDiscrete>`).

There are also methods for controlling addition of new individuals
within the species; `germination_probability` computes the probability
of germination, given an environment, `add_seeds` adds a number of new
individuals as seeds.

## Patch

A patch contains multiple species.  Direct interactions only happen at
the level of the Patch, mediated by the `Environment`, which `Patch`
controlls.  `Patch`es are far more autonomous than `Species` are.

In the deterministic model, there is just one idealised patch.  In the
stochastic model, there are potentially a large number of patches.

The method `size` returns the number of `Species` within a `Patch`.

It also looks after working out what the `Environment` (see below) is
like.  This is triggered:
  - on initialisation via `initialise`
  - when ode values are set with `set_ode_values`, as this changes the
    size of plants in the population.
  - after adding seeds via `r_add_seeds`
  - manually via `r_compute_environment` (NOTE: why is this needed?)
  
Adding individuals is a bit over-complicated at the moment.

  - `add_seeds` takes a vector of seeds (the same length as the number
    of species, given by `size`) stochastically germinates the seeds
    (that is, checks to see if they will survive in the current light
    environment) and then adds them to the appropriate species.
  - `add_seedlings` takes a vector the same length as the number of
    species and adds new individuals *without* checking if the seeds
    will germinate.  This is primarily used for set-up/testing.
  - `add_seedling` takes an index corresponding to a species (i.e., on
    [0,`size()-1`]) and adds a single seedling, without checking
    germination. This is used for the EBT where individuals do not
    stochastically die, but instead we track the probability that they
    did die there as part of their lifetime survival.

Births and deaths (in the stochastic version only) are handled via
`births`, which returns a vector of numbers of seeds produced, and
`deaths`, which loops over all species/individuals and checks for
deaths (NOTE: in contrast with `Species::deaths` we ignore the return
value -- should `Species::deaths` return `void` instead?).

Because Patches are somewhat self-sufficient, there are some functions
`step`, `step_deterministic`, `step_stochastic` for running the
system.  However, I might move these out into a controlling class
because they do complicate things.

## Metacommunity

A `Metacommunity` is composed of multiple `Patch`es.  At the moment
these are assumed equally connected (i.e., island model for
dispersal).  All seeds are immediately contributed to the global seed
pool (i.e., there is no "local" vs "global" dispersal distinction).

This class contains the main simulation control (via `step`, which
advances the simulation through whatever size step the ODE stepper
feels is appropriate).  It exposes (for now, at least) lower level
things, such as the `births` and `deaths` and `add_seeds`, all of
which behave like functions of the same name in `Patch`, except for
`add_seeds`, which also involves a *dispersal* step, where the given
vector (along species types) is dispersed among the `Patch`es.

## EBT

The `EBT` is the controlling class for the deterministic model (EBT
stands for "Elevator Boxcar Train", which is the numerical approach
that it implements -- perhaps this will change to a more informative
name).

It contains only a single `Patch`, and rather than stochastic arrivals
of seeds, there is a "schedule" of when cohorts are introduced
(`CohortSchedule`).

The EBT also needs the concept of "seed rain" (via the `SeedRain`
class), which keeps track of the total input of seeds.  This may find
its way into the stochastic version too (issue #28).

### CohortSchedule

This is effectively a queue of introduction times and the index of the
species being introduced.  The EBT just reads off the time that the
next introduction happens, advances to that time and introduces a new
cohort.

### Seed rain

The seed rain is simply the constant rate of seed influx.  This may
vary by species, so the `Patch` needs to make sure that the
`Environment` returns the correct seed rain for the current species
(see issue #22).  This is all a bit of a hack, and something that
could be fixed up, I think.  The seed rain is central to how the EBT
runs -- importantly changing the seed rain changes the initial
conditions.  As such, we can't change the seed rain of a running
population.

## Other biological components

### Environment

The `Environment` class pulls together everything that individuals
need to know about their environment.  It includes the canopy openness
at a given height (via `canopy_openness`), patch survival (via
`patch_survival`, but still flakey?), and the rate of seed rain (via
`seed_rain_rate()`).  The last two are only relevant to the
deterministic version at this point.

### Disturbance

Patches have a disturbance regime.  The model that the original
(Falster) version uses can be solved analytically, so we can tell when
a patch is created at what point it will suffer disturbance.

### Parameters

There are a handful of model parameters that don't apply to individual
`Plant`s (c.f. `Strategy`) -- these are combined into `Parameters`.
This includes a vector of resident plant strategies, but also
information on patch size, rate of light extinction, dispersal, and
disturbance.  `Parameters` also contains a number of control
parameters (that affect numerical routines, but not biological
processes) in the class `Control`.  This class is otherwise very
basic, and exists only to hold these pieces of data.

### Control

This is another data-holding class that contains parameters that
affect the numerical algorothms (e.g, the desired accuracy of the ODE
stepper, the quadrature routines, etc).  All the data members are
public.

# Utility components

## Lookup

A few classes inherit from the class `Lookup` (`Parameters`,
`Strategy`, `Control`, `ode::OdeControl`).  This exists purely so that
we can easily modify things from R.  So, `Parameters` contains many
(about 30) parameters that are all `double`, and that we will want to
be able to query and set from R but I don't want to make 30 different
pairs of get/set functions.  At the same time, these are not meant to
be modifable by any other classes and (in the case of `Strategy` at
least) private.

This defines a method `get_parameters` that returns a named R list of
all parameters within the object.  It also defines the method
`set_parameters` that takes a list containing new values for a subset
of parameters and modifies them within the object.

Classes that inherit from `Lookup` need to define a method
`do_build_lookup` that constructs the lookup table (constructed as a
`std::map<std::string, double*>`) and *may* define a method
`set_parameters_post_hook` that will be run after parameters are set.
This method can be used to compute constants/derived parameters.

## PtrWrapper

Creating a plant needs a *pointer* to `Strategy` object (this is
necessary as the `Strategy` outlives the `Plant`, and we don't want to
store 10,000's of `Strategy` objects).  Because of the
design/limitations of the Rcpp modules interface, we can't pass in a
parameter object to take the address from (even if we could, this
would do badly because it would be copied along the way and then the
address is invalid, or it might get destroyed before the plant).
Similar problems will affect the `Population`, as the `Metapopulation`
generally has control over this.

To get around that, I'm storing the `Strategy` within a small wrapper
class `PtrWrapper`.  This almost certainly duplicates `auto_ptr` or
something, and I should look into that (though our reference counting
semantics are a bit simpler).

## Splines

## ODE solver

We need an ODE solver that can cope with changing dimensions, as both
the individual-based and deterministic solvers have rapidly changing
dimensions (and the GSL solvers will have to reallocate *every* step,
which will become annoying).

The current implementation is a rewrite of the GSL solver, but using
vectors for allocation, which allows fairly graceful resizing.  The
stepping algorithm (Runge-Kutta Cash/Karp) does not depend on
long-term observation of the variables, so resizing space is all we
need to do.

### ODE State

The other major complication is that in contrast with the usual
presentation ODE problems (a pure derivative function doing (t,y) ->
(dy/dt)), we have a model that contains state (y) and can compute
derivatives (dy/dt).  So there ends up being two places that the state
and rates need to be.

The basic pattern is:
1. The ODE stepper asks the model for its current state `y`.
2. The stepper then advances the model by computing `dydt` for various
new times `t` and states `y`.

At the moment, this interaction is implemented using iterators for `y`
and `dydt`.  This is so that basically the same interface can be used
for each sub-component of the model; so long as the overal system
knows how long it is (dimension of `y`) we can just feed iterators in
that will distribute the values over the sub-components appropriately.


### Splines

We use splines for the light environment, based on the GSL library.
I've built a very simple class wrapper around the
construction/destruction/evaluation of the GSL splines.

There is a routine to adaptively create splines to resolve the light
environment to a desired accuracy (`AdaptiveSpline`).  

### Integration (quadrature)

Gaussian quadrature is used to compute the photosynthetic gain over a
canopy.  I've used the quadrature routines in the GSL library,
wrapping construction/destruction and evaluation.

### Root finder

This is needed only for the initial seed size calculation, which
itself is carried out only at the beginning of a simulation.  So
doesn't need to be anything too crazy.  It's also a wrapped version of
the GSL routines.

# Coding conventions

In most places I've tried to avoid abbreviation of variable names.

One place where short variable names persist is `Strategy`, where
variable names basically mimic those from the paper.  The scope of
these names is basically restricted to `Strategy` (which does nothing
with them) and `Plant` (which does a lot of calculations based on
them).

Many classes naturally contain others (`Species` contain `Plant`s,
`Patch`es contain `Species`, etc).  For these classes, I've used
`size` for the size of the container (so `Species`' size should
indicate how many individuals are there).  These all return `size_t`.
However, note that Rcpp translates that to numeric, and not integer
(to get the 64 bit precision, I think).

I've got `r_at` methods for accessing individual items from objects
that contain them (`Parameters -> Strategy`, `Species -> Individual`
and `Patch -> Species`, etc).  In the case of the Individual-related
bits, these can't go into the base classes because they require
templating.  From the R side, they use 1-based indexing and this is
bounds checked.

There are a number of `r_` wrapper functions that just bypass the
public/private interface for testing.  Options here could be to either
make these `protected` and inherit, or declare these inline rather
than having four lines where one would do:
```
ret_type r_method() const { return method(); }
```
vs.
```
ret_type r_method() const;
...
ret_type Class::r_method() const {
  return method();
}
```

The idea will be that a `r_` function should never be called except
for in a `r_` function or in `interface.cpp`.  All `r_` functions must
be public, but may just be direct wrappers around private functions.
If a public function has exactly the required interface then there is
no `r_` function wrapper.

In addition to testing, R-specific methods may be needed to check
arguments in a way that the compiled code need not do (e.g., checking
bounds before accessing containers), or for interacting with the
random number generator.

# Other notes

## Mutability

I've generally decided that parameters should not be mutable.  That
stops a huge number of corner cases, and should be a reasonable thing
to do as construction of the objects themselves is cheap relative to
how often parameters will be changed during a simulation.  So we just
start a population off with some types and run them.  If we want to
make a population with more types, we need to make a new population.
