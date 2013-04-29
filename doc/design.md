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
 
# Shared components

## Individuals

The "individuals" in the individual based model are very similar to
cohorts in the deterministic model; they have a state (mass, leaf
area, etc) and are capable of computing their physiological rates
given a light environment.

The key things that an indiviual needs to be able to:

* `get_height`, `get_mass_leaf`, `set_mass_leaf`: return its height and
  leaf mass, also set the leaf mass. [1]
* `compute_vars_phys`: given a light environment, update a number of
  physiological variables.  Could be seen as
  `update_with_environment`.
* `offspring`: return the number of offspring the plant is able to
  produce right now.
* `died`: test of the plant died right now.
* `germination_probability`: If the plant were a seed (assumes that
  the mass has not been changed)
---
[1] height is only used to determine the upper limit of the light
environment that we would care about; declarative programming would do
away with that.  The leaf mass could as easily be written as "plant
size" as nothing outside of the plant cares that this is the leaf mass
and not something else.  In fact, it might be nice to combine the two
and have height be the only interacting variable (check to see how to
derive leaf mass given height).

In addition, we need some methods by virtue of being able to compute
derivatives using `ode::OdeTarget` (see below).

* `ode_size`: return number of underlying variables
* `ode_values_set`: set underlying variables
* `ode_values`: return underlying variables
* `ode_rates`: return rates of change of underlying variables

The typical call by a controlling ODE stepper will look like:

```c++
ode_values_set();
env = compute_light_environment(); // given the new sizes of plants, etc
compute_vars_phys(env); // compute rates given light environment
ode_rates(); // get the rates of change
```

The most basic individual is `Plant`, which is a simple-minded plant.

In addition, there are some specialised versions of `Plant` that
should behave *(essentially) identically* when used in aggregate.
These are:

`CohortDiscrete`: if multiple individuals are born at the same time,
then we only have to compute their properties once (growth is
deterministic).  This simply adds a constructor that takes an integer
argument indicating how many individals are in a cohort, and modifies
a few methods where the output depends on the number of individuals in
the cohort.

`PlantApprox`: this version does the physiological calculations via a
different backend that uses an spline over the size axis (using
`PlantSpline`).  It changes very few methods, too.

`CohortDiscreteApprox`: this is a combination of the two approaches
(not yet implemented).

There will be a related individual type `Cohort` that is the
individual type for the EBT approach.  It will probably be generated
via composition though, with two individuals (top and mean).

#### Internals

There are two internal pieces that belong to `Plant`:

An internal class that contains the size variables.  Nothing else
knows about this, and it is defined within the `Plant` class
entirely.

A `Strategy` class that contains parameters that the `Plant` uses.
Other pieces of the program will push these around, but nothing is
allowed to see inside of `Strategy` except `Plant` (everything except
constructors are private, but the class has a `friend` relationship to
`Plant`.

#### More on Cohorts

More information on the `Cohort`s for the EBT:

Previously, the key bits of cohort change were:

* Mass of average and top individual (mu, bound)
* Density of individuals at the first time
* Survival at the top of the cohort
* Lifetime seed production at cohort boundary
* Seed output
* Number of individuals
* Density of largest individuals at the first time

We must be able to push new cohorts onto a list.  Probably to match
the individual based model, we should push on the back.

#### More on PlantApprox

The deterministic part can be done more carefully; suppose that we
have a vector of sizes of plants (at the patch level).  We can then
construct a spline over this size so that we can evaluate the target
functions (e.g., growth rate) as a function of size given the current
light environment.  Because interpolation will work well for
individuals that are close in size, this will greatly reduce the
amount of data we keep around, and how much work we have to do.

Most of this can be extended to the actual ODE calculation; this is
all deterministic too, so these just get evaluated on the splines.

Then we have some big lists/matrices of variables 
(list< vector<double, 3> > probably), and update those for every
plant, as before.

This has an advantage over the CohortDiscrete approach in that it
could allow for a true continuous time more easily (still would in
theory require that every birth and death reset the time, but I think
it could allow more flexibility there).

To do this, we would need an entirely new `Species` class; it would
have data members

```
class Species2 {
  std::list<Plant> representative_plants;
  std::list<double> mass_leaf; // main size axis
  std::list< std::vector<double> > other_variables; // mortality, fecundity
};
```

Working out how and when to refine / coarsen the grid would be hard.
Probably for large individuals we just want to keep all individuals
(so that the representative_plants would actually correspond exactly
to plants at some point).  But without rebuilding the spline every
time it will be hard to know.

Triggering `Plant` injection into `representative_plants` would
require quite a bit of work reimplementing logic in `AdaptiveSpline`.
An alternative would be to define a size "sapling" where the
self-thinning is less extreme and interpolate the values below this
over a fixed size number of points.


#### Stochastic equivalence of Plant and CohortDiscrete

A population based on `Plant`s and one based on `CohortDiscrete` will
stochastically diverge because of the differences in accounting in
`Species::died`; a `Species<Plant>` must go through and do independent
Bernoulli trials whereas a `Species<CohortDiscrete>` can do a binomial
draw when there is more than one individual in a cohort.

We could have a "hard mode" that forces Bernoilli trials; it would
look like

```
  bool CohortDiscrete::died_hard_mode() {
    const double m = mortality();
    const int n = n_individuals;
    for ( int i = 0; < n; i++ )
      if ( unif_rand() < m )
        n_individuals--;
    return n_individuals == 0;
  }
```

With this, the two cases would be stochastically identical.
However, at some point the ODE solver will do something with
scaling the error estimates by the number of variables (not sure if
it does so already, but if it's done on an average that will
happen) at which point we'll lose equivalence.

Another downside of this is that this would require more parameters
set somewhere.


## Species

Above the level of species, there are two places that "species"
appears as a concept:

1. as a structure/class of parameters, corresponding to a different
   strategy.
2. as a book-keeping device, collecting individuals within a
   population
   
## Patch

A patch contains multiple species.  In the deterministic model, there
will be just one patch.  Previously this was population.

It's not totally clear how seed addition should look here:

* The `add_seeds` function will always test to see if seeds have
  germinated and survived dispersal.

* From R, sometimes we want to just dump seeds in and have them skip
  germination entirely.  In this case, we work via the
  `r_add_seedlings` function, which still takes just an integer vector
  and still adds individuals in as if they were seeds.  However, every
  seed germinates with probability 1.
  
Note that this distinction is made only in `Patch` and
`Metacommunity`, because `Patch` is where germination probability is
computed and germination tested.  For Species, there is no
`add_seedlings`, and `add_seeds` just adds seeds.  It could be updated
for consistency with `Patch`, but I think `add_seeds` is possibly
clearer.

### Disturbance

Patches have a disturbance regime.  The model that the original (EBT)
version uses can be solved analytically, so we can tell when a patch
is created at what point it will suffer disturbance.

When should this happen?  Perhaps we should do `step_deterministic` /
check if we've passed death time / `step_stochastic`.  This way we
will definitely complete a step (deterministic) but not produce
seeds.  This approach very slightly increases the rate of
disturbance.

A more exact approach would be to have the system "know" about what
times are special.  Then we just refuse to run past that time.  This
would require another set of bookkeeping at the level of
Metapopulation, but we will need this (somewhere) for the EBT version,
too.

Then, when we step the Patch, we should check about disturbance.  On
the other hand, we might want the Metacommunity to decide when the
patches get disturbed, so I'm not totally sure about the long term
plan here.

## Metacommunity

A metacommunity is comoposed of multiple patches.

At the moment these are assumed equally connected (island model for
dispersal).  This extends to seed producing patches, with all seed
immediately joining a global seed pool.

## Parameters

We need a pointer to all the parameters that we can use to seed new
populations (containing things like a vector of all species
strategies), and a few parameters that apply at the level above
species (light extinction, disturbance regime, etc).

The nice thing about this approach is that we can centralise ownership
and lifespan of the parameter pointers.

Need to decide on names for accessor functions.  Probably change
`set_params` to `set_parameters` to continue avoiding abbreviations.

## Standalone Plants, Populations, etc.

Creating a plant needs a *pointer* to `Strategy` object (this is
necessary as the `Strategy` outlives the `Plant`, and we don't want to
store 10s of thousands of `Strategy` objects.  Because of the
design/limitations of the Rcpp modules interface, we can't pass in a
parameter object to take the address from (even if we could, this
would do badly because it would be copied along the way and then the
address is invalid, or it might get destroyed before the plant).
Similar problems will affect the `Population`, as the `Metapopulation`
generally has control over this. 

To get around that, I'm storing the `Strategy` within a small wrapper
class.  This almost certainly duplicates auto_ptr or something, and I
should look into that.

It's possible (probable?) that for `Species` and above there is
absolutely no benefit from this distinction, and dealing in plain
objects will be better.

Also mention when doing the cross of PlantApprox and CohortDiscrete
that I'm trying to avoid DDD.  Traits might be better?

So, define a class member `standalone`, which will be false normally
but true if a `Plant` is created without a `Strategy` pointer.  In
that case, it will either create the default `Strategy` or take a copy
of a provided `Strategy` object, and be responsible for cleaning up on
deletion.  Copies involving a standalone `Plant` will *copy* the
`Strategy` object and create a new "standalone" `Plant`.  This is done
through the helper class `WithStrategy`.

A slightly better approach would be to template a simple smart pointer
that did all this, perhaps?  Same as existing case, but with
`Strategy` as a template parameter, and save the pointer in field
`ptr`, then on use do

```
class Plant : protected SimpleWrapper<Strategy> {
  Plant(Strategy s) 
    : SimpleWrapper(s),
	  strategy(ptr) {}
};
```

but that leaves us in the same position with copying as before I think
(i.e., strategy will point at the wrong pointer).  But if we template
the base class and have concrete uses of them, then this goes away.

## Dispersal

I might try and do dispersal through a map object with the pointer to
the relevant strategy and the number of seeds.  Alternatively, we
could just use the index, as that is known at the correct level.

So, the metapopulation will create a "seed pool", which is simply an
integer of the appropriate length.  It then disperses the seeds.

Where does the probability of dispersal belong?  Is this a property of
a plant, or is it a property of a patch?  It seems that there is a
possibility for different dispersal probabilities to occur at
different points; one leaving the plant, one in the source patch, one
in the destination patch, and one depending on the distance between,
etc.  For now, I'm moving it into the patch.  I'll decide if it's paid
by the source or destination patch once we get a metapopulation
going.

Eventually, I might need a concept of "incoming seed" within a Patch;
this would allow us to increment this during the dispersal step.  This
means that a vector of "types" is probably an appropriate intermediate
data structure to go through.

One way of thinking about where the "correct" place to compute things
is what happens when we go to a spatial model?  The patch will say
"what is your germination probability, little seed" for each seed
(within each species) and append it as needed; it really won't happen
at the species level.

## Seed rain

In the EBT version of the model there is the idea of "seed rain",
which I have yet to add in.  Need to think about how this is done.

## Deaths

## Births & Deaths

Not sold on the names for the birth/death functions.  They are all
going to return different types and work in different ways (some
return bools, others modify lists and return integers, etc).

## Adding strategies, constants and standalone issues

Where should strategies be added?  The two most basic options seem to
be:

1. In population
2. In the parameters only

At some point, constants do need sorting out.  So, perhaps on
`add_strategy`?  I.e., we do `add_strategy(x)` (or `set_strategy`),
and then make a standalone plant which computes the constants?  Should
be a better way of doing this.  Perhaps a static method for Plant?

### Changing plant parameters

Probably don't allow plant parameters to be modified?  Otherwise they
change the paramters of all individuals, which will rarely be what is
wanted.  Exception for standalone plants?

### Exception safe code

There is a good chance that the root finding can fail, but it is run
on construction (e.g., standalone Plant).  Probably not a good idea.

## Constants

There are probably still some constants that are settable at the level
of `Strategy`.

For example, `(c_acc * s)` appear together and only together in
computing `fecundity_rate`.  As such they are entirely colinear.

In computing `mortality`, the value
`strategy->c_d0 * exp(-strategy->c_d1*strategy->rho)`
is constant.

There are others in `compute_vars_size`:
  - `rho / theta * a1 * eta_c` (`sapwood`)
  - `rho         * a2 * eta_c` (`heartwood`)
  
When and where do we compute these?  At the moment they are set in
`Strategy`, but this splits logic over too many files.  The advantage
of this approach is that we could set it up to run after every
parameter set.  The disadvantage is that to do anything complicated
(such as the seed leaf mass calculation) we will have a recursive
reference.

Probably better to move the constant calculation within the `Plant`
class.  It might be desirable to have this be a static method?

So then, how do the constants get computed?  Who creates the
Parameters instance?  In the full model it will be done by
`Parameters` generally, but for standalone cases the `Plant` itself
takes care of this.

### Accessing constants

Currently there is no way of accessing the constants (especially
initial seed mass).  How should we do this?  Should we do this?

### Default values

I have now had to add default values for the four core traits, too,
otherwise we get cryptic errors when trying to compute the initial
leaf mass.  It doesn't really make much sense to make a strategy that
doesn't know it's lma, so this seems OK.

### Initialising constants

This should be done by anything that takes a `Strategy` object (rather
than a pointer).  It will be copied, so there is no danger.  This will
be the same as things that are "standalone".  However, when a
`Parameters` is made, can we assume that this has taken care of
sorting this out for us?  In general, this issue is a bit of a design
wart.  But given the dependencies and trying to avoid a circular
dependency.  Not sure how I can do this better though.

## Functions within plant

I've pulled into separate functions the physiological variables that
depend *only* on size variables (for now at least).  The reason for
this is that things like `net_production` depends on `assimilation`
and `respiration`, and only makes sense to be called after these.  So
by the time it is callable we already have the values stored in the
class variables.  Ditto for `growth_rate`, which depends on
`leaf_fraction`, `net_production`, `leaf_fraction`.

In some ways, a slightly nicer way might be to push these into a
separate method `compute_allometry()`, but some of these are only
worth computing if net production is positive.

`assimilation`/`compute_assimilation` is the odd one out here, as it
depends on the environment variable.  It may be the case that this
ends up set for the individual (on entry to `compute_vars_phys`,
perhaps), in which case this difference disappears.  But that would
really require that it is used in something other than computing
`assimilation` to be worth the confusion.

For the same reason, the individual allometry functions are not
written separately; each part depends on previously calculated
things.  So this is all done in `update_vars_size()`.

## Utilities

### General

While I now basically have all of the bits that we need.  However,
there is a diversity of interface issues still remaining that could be
tightened:

Integrator and FindRoot take their (modest list of) control parameters
during initialisation.

ode::Solver has a whole separate class for control parameters (not all
of which appear to be enforced yet or in the right place --
e.g. `step_size_min`, `step_size_max`, and `no_steps_max`).  However,
there is no ability to set them to anything other than the defaults.
There are rather a lot of possible parameters too.

`Spline` doesn't actually have *any* control parameters, apparently.
However, `AdaptiveSpline` does, and that needs harmonising.  It also
has a modest list, so can be done via the constructor.

Bounds setting: `AdaptiveSpline` sets that bounds and target
separately (and uses the old target style).  Move this over to the
same approach as the other utilities.

One way of doing this is to derive from a `WithTolerance` class that
takes care of getters, setters, and defaults?  Could just be more
abstraction than is helpful though.

On the other hand, constructor arguments are often shunned, so perhaps
all these bits should take an `init()` method or something?  The
`set_control` method to `AdaptiveSpline` but with error checking is
possibly the trick.

The Integrator does the integration and the bounds all at once.  That
is probably the way that I should use the adaptive spline too.

### Organisation

Not sure where the utilities code should all go, from a namespace
point of view.  I have `Spline` and `AdaptiveSpline` both within
`spline`, `Solver` is in `ode` (along with its support)

`solver_ode.{cpp,h}` filename is inconsistent with everything else.

### ODE solver

We need an ODE solver that can cope with changing dimensions, as both
the individual-based and deterministic solvers have rapidly changing
dimensions (and the GSL solvers will have to reallocate every step,
which will become annoying).

At the moment, this is implemented using iterators for `y` and `dydt`.
This is so that basically the same interface can be used for each
sub-component of the model.  Alternatively, we could work with (say)
`Patch::derivs` taking `y` and `dydt` as `const` and non-`const`
references respectively, and then do the rest of the inner workings
with iterators.

However, I do need a typedef somewhere in the ODE so that I can do a
`ode_iter` and `ode_const_iter` and skip some of the ugliness that I
currently have.

There are a lot of duplicated functions (I see `ode_values`,
`ode_rates`, `ode_values_set`, `derivs` *and* the R versions of these.
If we inherit from an interface class we could simplify a lot of this?

Note that the `derivs` one might be hard.  In particular, for `Plant`,
we need the light environment to compjute the derivatives, and that is
not stored as part of the class.  So either we modify `Plant` to take
the light environment (which simplifies the `FunctorBind1` hack a bit)
or we just can't do this.  The downside of this is that the logic of
when to compute physiological variables becomes hard.  Is it on set?
Is it on request?  What if the the value is NULL/expired/etc.

#### ODE State

The more clever ode systems do things like tell you if a solver is
capable of making use of the exit values from one cycle to be the
input values for the next.  I'm trying to avoid doing too much of this
sort of optimisation right now because I know that the problem
dimension will change pretty much every point.

The problem is that we must easily update the state; we're not
generally allowed direct access to the ODE state object, but perhaps
we can arrange to set it with an iterator? 

What I'm trying to avoid is a situation where we go through and build
a new vector just to get destroyed being passed into `set_state` (this
is what is currently happening).  However, it's not really clear that
there is a better solution.  At the moment, the function used is an R
interface function from ode::OdeTarget, which is really not great.

In contrast with the other copies, this is going to create a large
heap-allocated copy and then discard it straight away.  Hopefully the
compiler is smart enough to deal with that.  It's possible that we
could use an intermediate piece of storage at the level of the class.

Simply having `ode_state` take an iterator would solve at least one
copy, but not be hugely useful as the we're still creating a temporary
object.  If the ode system provided an iterator to its internal
storage we could fill that, but that's ugly too.

### Splines

We use splines for the light environment, but I have some code I'm
quite happy with that works here.  However, the reallocation is slow,
so this is likely to find itself ported over from GSL to proper C++ at
some point if it persists in being slow.  However, we should be able
to do that entirely transparently.

TODO: Crash if eval is called before init.

### Numerical integration (quadrature)

This has ended up being a very simple interface, and I think that I'll
use this as a model for the other utilities where possible.

### Root finder

This is needed only for the initial seed size calculation, which
itself is carried out only at the beginning of a simulation.  So
doesn't need to be anything too crazy.

# Compilation issues

Rcpp makes compilation very painful, so we really want to load that as
little as possible.  However, the utilities depend on it, and whether
it is included causes Rf_error() vs. error() differences (amongst
others).

# Filenames and conventions

Deciding if to go with class variables with underscore (google style
guide).  The large number of class variables and equations will make
this a pretty pervasive change.

Initially, I used single letter variable names in `Plant` -- for
consistency with the paper.  This is complete for the variable names,
but functions still have a way to go.

One place where short variable names persist is `Strategy`, where
variable names basically mimic those from the paper.  The scope of
these names is basically restricted to `Strategy` (which does nothing
with them) and `Plant` (which does a lot of calculations based on
them).  So I should set all the data content of `Strategy` private and
declare `Plant` to be a friend.

The downside of the friend approach is that there are things that I
really don't need `Plant` to be able to see (`lookup_type`,
`build_lookup`, `lookup_table` and `lookup`).  All of these are based
around the issues of the data lookup.  So provide a virtual base
class, privately inherit from this, and then the issue goes away?
That has the nice effect that the lookup problems go away.

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


# Consistency

init vs contructors
reset, clear, etc.
helper vs wrapper

Need to decide how the ODE-related stuff gets named.  I like `size`
for the size of the container (so Species' size should indicate how
many individuals are there).  But at some point we might have to
change that?  I think that the ODE system just uses the size of the
created state, and we can get that with `push_back` if needed?

Note that all `size` should return `size_t`.  However, note that Rcpp
translates that to numeric, and not integer (to get the 64 bit
precision, I think).

I've got `r_at` methods for accessing individual items from objects
that contain them (Parameters -> Strategy, Species -> Individual and
Patch -> Species, etc).  In the case of the Individual-related bits,
these can't go into the base classes because they require templating.

It's possible that we could inherit from a container class that wraps
stuff up, but that's a lot of mucking about for a small amount of
duplicated code (plus it won't work for species, as that uses a list
not vector for storage anyway).

## Public, private & R interface

I want to test more functions than I actually want to make public, I
think.  For example, `Patch::height_max` is not exactly bad to expose,
but it's not really needed to define the interface either (it's
essentially an implementation detail that only height-related models
care about).  However, I do use this a bunch from R.  Perhaps make it
private and provide an `r_height` method that I can access from R to
make this separation clear?

The idea will be that a `r_` function should never be called except
for in a `r_` function or in `interface.cpp`.  All `r_` functions must
be public, but may just be direct wrappers around private functions.
If a public function has exactly the required interface then there is
no `r_` function wrapper.

The ordering within the header file remains public/private, but within
the cpp file, the public `r_` functions are all kept at the bottom of
the file.

One reason for this is that many input-taking R-exposed function needs
to check arguments in a way that the compiled code need not do (e.g.,
checking bounds before accessing containers.

# 0-based vs 1-based indexing

There are a few R-based functions that modify things.  At the moment
they all use C-style 0-based indexing (with bounds checking).  They
should probably all move over to useing R-style 1-based indexing.

# Callbacks and functors

There are two basic ways that callbacks are being done --
AdaptiveSpline does it in an ad-hoc way, but the root finding (and the
quadrature approach) do the call backs with functors (I guess that the
ODE class is doing it a third way!).

Move everything that make sense over to functors.

## Update:

I apparently did not know about `std::bind*` and `std::mem_fun`, which
will do a lot of this for us.  This will probably make the code a lot
less fragile.  However, it does all work at the moment, and will take
some care to fix.

# Alternative models

The Kohuama and Takada model would fit quite easily into this
framework, and might be a nice way of testing some ideas about the
"right" level of generality.  It would generate a different looking
"environment" object (sums of leaf area over discrete strata), but the
idea for plants would still be to convert this environment into
growth.

The "size" value becomes the partition that a species is in for the
individual based model.  For the EBT-style model we are also tracking
the fraction of the population in this size class (but unlike the EBT
we never bother splitting cohorts).

# Compilation issues

In theory, we can trim out the Rcpp bits with some ifdef'ing.  Will
require some additional includes and coping with bailing.  Would help
enable model to be embedded in other control forms (e.g. python with
boost.python).

To that end, all Rcpp-involving functions should be r_ prefixed.

With clang 3.1 (osx 10.8) this is my system Makevars
```
CXX=clang++
CC=clang
SHLIB_CXXLD=clang++
CPLUS_INCLUDE_PATH=/Library/Frameworks/R.framework/Versions/2.15/Resources/library/Rcpp/include/
CXXFLAGS=-O2 -Wall -pedantic -Wconversion 
```

There are more `(int)` conversions that desirable, mostly working
around Rcpp at the moment.
