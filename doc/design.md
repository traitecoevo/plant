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

## Individuals (Plant, etc)

The "individuals" in the individual based model are very similar to
cohorts in the deterministic model; they have a state (mass, leaf
area, etc) and are capable of computing their physiological rates
given a light environment.

The key things that an indiviual needs to be able to:

* `get_height`, `mass_leaf`, `set_mass_leaf`: return its height and
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

### Key variables

There are a few key variables -- these are plant size (leaf mass),
mortality and fecundity.  To make the ODE stepping work, we need to be
able to get and set these variables as well as get the rate of change
of these.

We also need height in a couple of places.  Given that there is a 1:1
correspondance between height and seed mass, it might be worth
reworking all the calculations that way.  The only place that this is
actually hard is getting dh/dt from d(mass_leaf)/dt.

#### Abstracting these

There is quite a bit of duplication floating around with the key
variables.  Every one has a get/set/rates function, plus two spots
somewhere where numbers are stored.  It seems that they could be
paired together into some generic "ODE variable" type.

These variables are actually kind of weird; we'll often want to access
them from R (especially for tests), but never from patch etc.  So
something like

```
plant$variable("mass_leaf")
plant$variable("mass_leaf") <- x
plant$rate("mass_leaf")
```

could be a nicer way of capturing the dynamics and extend more easily
to different numbers of traits.

### Changing plant parameters

Probably don't allow plant parameters to be modified?  Otherwise they
change the paramters of all individuals, which will rarely be what is
wanted.  Possibly a change in parameters could trigger a clear of all
individuals?

### Functions within Plant

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

### Internals

There are two internal pieces that belong to `Plant`:

An internal class that contains the size variables.  Nothing else
knows about this, and it is defined within the `Plant` class
entirely.

A `Strategy` class that contains parameters that the `Plant` uses.
Other pieces of the program will push these around, but nothing is
allowed to see inside of `Strategy` except `Plant` (everything except
constructors are private, but the class has a `friend` relationship to
`Plant`.

### Default values

All parameters have default values following those in the paper.  The
four "core" traits have default values added too; these came from the
orginal EBT implementation (rationale is: it doesn't really make much
sense to make a strategy that doesn't know its lma).

### Constants

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
  
Currently, these are done after setting parameters within `Strategy`,
using a static method from `Plant`.  This keeps all the logic within
`Plant`, and `Strategy` is simply a container of parameters.

Currently there is no way of accessing the constants from R (or
indeed, anything except for `Plant`).

### More on Cohorts

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

There are two different sorts of cohorts that were previously rolled
together:

* CohortMean -- based around the mean value of a potentially large
  number of individuals within a cohort.  This is the "traditional"
  EBT model, and requires some care at the bottom individual.
* CohortTop -- based around the size of the largest individual.

The first of these involves more approximations, and should be less
accurate than the second.

#### CohortMean:

For cohorts *except* the smallest (i = 2, ..., k)

Define the number of individuals, $\lambda_i$ as
$\integrate_{\Omega_i} n(x, m, a) \mathrm{d}m$.

The rate of change in the number of individuals is
$-\integrate_{\Omega_i} d(x, m, a) n(x, m, a) \mathrm{d}m$
which is approximately
$-d(x, \mu_i, E)\lambda_i$.

This is the only extra piece of work.  The fecundity rate is as
before, and survival is through the death rate, as before.

This probably wants to be a special case of Plant, as we'll just use
the machinery directly.  That must be easier than having to set up
forwarding methods for everything.

#### CohortTop

* Size increases simply with the growth rate.
* Mortality "variable" increases with the mortality rate (survival is
  exp(-mortality))
* Density of individuals (n(x,m,a)) requires an extra trick or two.
  $n(x,m,a) = n(x,m_0,a_0) exp(-\integrate_{a_0}^a X(a') da')$
  where
  $X(a') = \partial growth_rate/ \partial m + death_rate$ -- we can
  integrate that X(a') over time with the ODE stepper.
* Seed production: $pi_0 * pi_1(x, m_0, a_0) \integrate_{a_0}^a F(a') da'$
  $F(a') = fecundity_rate(a') * S_I(a') * S_P(a')$

To compute the partial differential equation, suppose we have a
function that returns growth rate for a given mass.  Then, we can bind
that with a Functor and pass it to a gradient finding function.  This
would be a target for eventual replacement with a better adaptive
gradient calculation function, or with calculation using the
approximate plant idea.

The survival calculation is a bit tricky.  In the document that Daniel
has made, we have see production has

$$
R(x, a_0, a) =
\pi_1(x, m_0, a_0)\int_{a_0}^a \pi_0 * f(m(a_0, a'), a') *
  S_I(a_0, a') * S_P(a_0, a') d a'
$$

where $\pi_1$ is the germination probability.  Note that

$$
S_i(a_0, a) = \exp(-\int_{a_0}^a d(m(a_0, a')) d a')
$$

which is simply the probability not dying between $a_0$ and $a$.  So
$1-\pi_1$ could become the initial condition for that.

In some ways it would be good to hijack the seed production bits
within Plant, but we don't because that could get a little confusing.
So we use a new variable `fecundity_survival_weighted`.

Computing the initial survival needs to be done as a cohort is added.
It is not necessary as it is *created* though, and initial survival is
set to 1 on creation.

There is a real ugliness in this at the moment: to compute the seed
output, we need to know the *patch* survival.  This is something that
we never need to know for the individual-based model.  This means that
`compute_vars_phys` needs its signature updated, but we can't overload
it because it's a virtual function.  At this point, I've called it
`compute_vars_phys_surv`, as it's really both physiology and
survival.  This will hopefully change.

##### Initial conditions

*This section should probably evolve with the documentation in EBT.md*

There are a series of cohorts -- these are at masses $m_1, m_2, ...,
m_n$ for $n$ cohorts.  However, to compute things over the
distribution of cohorts (such as light with equation eq:light), we
also need a "pretend" cohort $m_0$, which always has the seed mass
(i.e., it can never grow).  To compute the *density* of this cohort,
we will need to compute its growth rate: equation eq:boundN (&
following eqn).

At the moment, cohorts are set up with:
  * `mass_leaf` = leaf mass at birth, given growth model
  * `fecundity` = 0 (trivial)
  * `mortality` = $-log(\Pr({\rm germination}))$ (i.e., probability of
    dying will be 1 minus probability of germination.
  * `density`   = $y_x / g$, where $y_x$ is the seed rain (set by
    environment?) and $g$ is the growth rate of a seed in the current
    environment.  Or 0 if $g \le 0$.

So, who looks after the seed rain?  It never changes during a
simulation run, depends most strongly on the Species, but is something
that isn't really "Strategy".  It feels like it *should* go into the
Environment (partly because that's the other argument to
`CohortTop::initialise`) but it is not ever used by the stochastic
version (but then neither is things like `births()` and `deaths()`
that `CohortTop` inherits from `Plant`.

One way to think about this is *where would we set the seed rain*.
We're likely to do something like:

> "run this simulation with Parameters 'p' and seed rain 'x'"

So, it probably ends up in the Parameters vector at some point.  That
makes sense because we spoke about having a seed rain vector there
anyway.

**Except**, that we actually add strategies to the parameter vector by
pushing them back.  This *feels* more like an environmental thing.
So, have the environment containing a vector of seed rain; this is
much like the vectors that we were using for seed in/out for patches
(and perhaps we could repurpose this for same).

1. Probably subclass the core plant bits out from the dispersal stuff
at some point.  This is hard with `Species` needing to use them
though.

#### Implementation notes

At the moment, this is implemented by inheriting from Plant, as we
need most of the machinery there.  However, there are some aspects of
`Plant` that are tied up with the stochastic version (`offspring` and
`died` in particular) that aren't needed in the EBT version.  However,
`Species` will also have these, so need to be careful about how these
are exposed.

#### Bottom cohort

When dealing with integrating over the distribution of individuals, we
need to include a cohort that is exactly at the seed mass and does not
grow -- essentially we need a bottom cohort that does not change.
Ideally the plant does not need to know about this at all, so we
really want to do this at the level of the species.

This affects only `compute_vars_phys_surv`, really, and the
derivatives calculations.  We just don't want to 

The trick in all of this is that to compute the boundary condition
$n(x, m_0, a_0)$, we need to know $g(x, m_0, a_0)$.

So: during `ode_values_set` we want to set *value* of the bottom
cohort, but we don't know $g$ yet because that requires the
photosynthetic environment (which is not known).  So far OK -- it will
be set to the previous value.

Then we *do* compute $g$ and can update $n$.  The problem is that the
ODE stepper will do that for us again later (will set it back to the
previous value because)

So: the basic idea (within the `derivs` function)

(once before doing the derivs step a number of times)
* System takes values from the current species -- this will include
  the mass as the initial leaf mass, survival of zero, mortality of
  the germination probability in the *previous* light environment for
  the initial conditions for the ODE solver.

* ODE solver *sets* values within the system.
* Compute the light environment
* Get all cohorts to update their physiological rates (including the
  lowest cohort).
* Collect the rates of change for everything.  However, for the
  smallest cohort, set the rates of change to zero.
  
We can do this by specialising the `ode_rates` method, but there are a
few other bits that need doing too:

* On creation of a species, or on `clear()` we need to create the
  bottom cohort.  This needs to involve creation of an empty light
  environment and a bit of a hassle.
  
The big issue is that the cohorts need to be able to compute patch
survival since their birth to be able to compute their seed output.
The parameters for this are sorted out in the Parameters object that a
Cohort does not have access to.

So, what I'm planning on doing is having a functor or something that
we can pass the Cohort on initialisation.  Of course, it's not clear
how we can initialise that without a new method (that would require
changing the Species class a lot).

An easier way of getting something working in the short term would be
to put the survival function right into the plant.  However, even then
getting the actual parameter values added is tricky.

A longer term solution could see something like the Patch survival
coming along as part of the environment.  So, with that as a longer
term aim, I've just got something fairly quick and dirty working now.

### More on PlantApprox

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

### Stochastic equivalence of Plant and CohortDiscrete

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
   
Some of the bits in `Species` aren't going to be wanted for
`Species<CohortTop>` for the EBT.  In particular, 

* The `leaf_area_above` calculation needs to be done carefully (but
  see `CohortDiscrete` for how that can be done simply).
* The `compute_vars_phys` calculations needs to be done via
  `compute_vars_phys_surv`, passing in the patch survival.
* `births` and `deaths` are just not really possible, nor do they make
  sense -- this is probably the biggest ugliness.
* `add_seeds` may need some care because of issues around getting the
  entry level cohort done correctly.

It might be best if we make a special class of Species -- something
like:

```
class SpeciesCohort : public Species<CohortTop> {
  ...
private:
  // methods to disable go here.
};
```
   
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

### Patches for EBT/individual based

Probably we will have a merge/split with patches -- there is a lot of
shared code.

* both
  - `step()`
  - `initialise()` (if kept?)
  - `height_max()`
  - `canopy_openness()`
  - `compute_light_environment()`
  - `compute_vars_phys()`
* stochastic only
  - `step_deterministic()/step_stochastic()`
  - `births()/deaths()`
  - `add_seeds()/add_seedlings()` (possibly different interfaces?)
  - `germination()`
* ebt only
  - `set_seed_rain()` (possibly?)
  - any code related to cohort timing/next step, etc.
  - mutant rerunning code (possibly?)

The big hassle will be that there will end up with some duplication of
exports for the interface code.

This is complicated by the fact that we have inheritance and
templating going on, interacting with the desire to export these via
Rcpp (which need concrete instances).  Inheriting off teplated classes
is a bit tricky for me:

```
template <class Individual>
class PatchStochastic : public PatchBase<Individual> {
  ...
};
```

won't work --
  http://stackoverflow.com/questions/3277812/c-template-class-inheriting-another-template-class-with-a-template-specified-i
This will make the rest of the code really ugly.  The issue here is
apparently "The main reason C++ cannot assume anything here is that
the base template can be specialized for a type later"
  http://stackoverflow.com/questions/11405/gcc-problem-using-a-member-of-a-base-class-that-depends-on-a-template-argument

This all suggests some fairly bad design problem.  I suspect that the
class is trying to do too much.  For now, I might just have to
duplicate the code in Patch for EBT and work back to sort it out
later.  A bit of a wart, but I'm not seeing a really obvious solution
here.

I feel that there's a core of code in there that both bits need, and
then there is a way of actually interacting with the Species within
the patch that will vary depending on how we use it (as such,
inheritance is probably inappropriate because we'll treat the objects
totally differently).  There may be a compositional way of doing this,
but I don't know how to do this following an open/closed model
(especially without exposing much of what is going on in Species).
It's possible that we can do this with passing back iterators or something?

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

I'm not totally convinced that the dispersal function should belong in
the class (the only class data that it uses is the number of patches).
It seems that a more generic dispersal object could be useful here.
However, it's difficult to see how that would work without knowledge
of the underlying spatial arrangement, which would cause some leakage
between classes, probably.

## EBT

This is done quite differently than the Metacommunity, because there
is only one "patch".  More differently even, is that the schedule of
arrivals is known in advance.

I'm unsure of where to put this.  Two ways seem most obvious:

1. Each species could look after its own introduction times.  This
   requires an additional field in `Species`, which is undesirable as
   it's not a self-sufficient class anyway (it can't be "run" because
   it doesn't know how to compute the light environment).  This also
   requires that we cycle through the species at every step to see
   what the next stopping time is.
2. The EBT will look after the times.  I'm going to go with that.

### More on dispersal

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

### Seed rain

In the EBT version of the model there is the idea of "seed rain",
which I have yet to add in.  Need to think about how this is done, and
how necessary it is in the stochastic version.

It would be good if seed rain was activatable at both the the Patch
and the Metacommunity; that way once seed rain has stabilised at the
Metacommunity level, we can easily run individual patches with that
level of seed (simulating a new patch that has no effect on the other
patches).  This brings the individual based simulation more closely in
alignment with the EBT.

In theory this is quite simple -- if the seed rain is a constant
process, unaffected by the model itself, then it is representable by a
vector of rates (r1, r2, ..., rn) for n strategies; over a period of
time dt, the expected number of seeds that arrive is then Poisson
distributed with mean ri * dt for the ith species.  We can just take a
draw of that many seeds and add them to the seed input from dispersal.

## Parameters

We need a pointer to all the parameters that we can use to seed new
populations (containing things like a vector of all species
strategies), and a few parameters that apply at the level above
species (light extinction, disturbance regime, etc).

The nice thing about this approach is that we can centralise ownership
and lifespan of the parameter pointers.

Need to decide on names for accessor functions.  Probably change
`set_params` to `set_parameters` to continue avoiding abbreviations.

## Control

In contrast with Parameters, this is where I'm thinking we can dump
all the algorithmic control parameters.  The separation is fairly
arbitrary, of course.  Not sure how we'll assign these, but thinking
something like the way that strategy is done on pointers seems about
right.  We can just provide one on construction?  Then bundle it
together with Parameters so that only Species & Plant need to bother
about getting one specially.

# Design issues:

## Standalone Plants, Populations, etc.

Creating a plant needs a *pointer* to `Strategy` object (this is
necessary as the `Strategy` outlives the `Plant`, and we don't want to
store 10s of thousands of `Strategy` objects).  Because of the
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

### Exceptions

There is a good chance that the root finding can fail, but it is run
on construction (e.g., standalone Plant).  Probably not a good idea? 

Probably I should think about moving from `Rf_error` to proper C++
exceptions.  Apparently `Rf_error` jumps over C++ destructors, so will
cause leaks.

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

The actual biology is in the `model` namespace, and testing functions
are in `test` namespaces (within `util` or `model`).

### ODE solver

We need an ODE solver that can cope with changing dimensions, as both
the individual-based and deterministic solvers have rapidly changing
dimensions (and the GSL solvers will have to reallocate *every* step,
which will become annoying).

The current implementation is a rewrite of the GSL solver, but using
vectors for allocation, which should allow fairly graceful resizing.
The stepping algorithm (Runge-Kutta Cash/Karp) does not depend on
long-term observation of the variables, so resizing space is all we
need to do.

The other major complication is that in contrast with the usual
presentation ODE problems (a pure derivative function doing (t,y) ->
(dy/dt)), we have a model that contains state (y) and can compute
derivatives (dy/dt).  So there ends up being two places that the state
and rates need to be.

At the moment, this interaction is implemented using iterators for `y`
and `dydt`.  This is so that basically the same interface can be used
for each sub-component of the model; so long as the overal system
knows how long it is (dimension of y) we can just feed iterators in
that will distribute the values over the sub-components appropriately.

#### ODE State

The more clever ode systems do things like tell you if a solver is
capable of making use of the exit values from one cycle to be the
input values for the next.  I'm trying to avoid doing too much of this
sort of optimisation right now because I know that the problem
dimension will change pretty much every point (and individuals might
die in the middle, not just at one end).

The problem is that we must easily update the state; we're not
generally allowed direct access to the ODE state object, but perhaps
we can arrange to set it with an iterator?

What I'm trying to avoid is a situation where we go through and build
a new vector just to get destroyed being passed into `set_state` (this
is what is currently happening).  However, it's not really clear that
there is a better solution.  At the moment, the function used is an R
interface function from `ode::OdeTarget`, which is really not great.

In contrast with the other copies, this is going to create a large
heap-allocated copy and then discard it straight away.  Hopefully the
compiler is smart enough to deal with that.  It's possible that we
could use an intermediate piece of storage at the level of the class.

### Splines

We use splines for the light environment, but I have some code I'm
quite happy with that works here.  However, the reallocation is slow,
so this is likely to find itself ported over from GSL to proper C++ at
some point if it persists in being slow.  However, we should be able
to do that entirely transparently.

### Numerical integration (quadrature)

This has ended up being a very simple interface to GSL routines, and I
think that I'll use this as a model for the other utilities where
possible.

### Root finder

This is needed only for the initial seed size calculation, which
itself is carried out only at the beginning of a simulation.  So
doesn't need to be anything too crazy.

# Filenames and conventions

## Variable naming

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

## Method consistency

Unresolved issues around:

- init vs contructors
- reset, clear, etc.
- helper vs wrapper

## Container-like classes

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

## Public, private & R interface

I want to test more functions than I actually want to make public, I
think.  For example, `Patch::height_max` is not exactly bad to expose,
but it's not really needed to define the interface either (it's
essentially an implementation detail that only height-related models
care about).  However, I do use this a bunch from R.  Perhaps make it
private and provide an `r_height_max` method that I can access from R
to make this separation clear?

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

### An alternative approach

Some of the classes get quite large because of all the R interface
code bolted on (e.g. `Patch`, in particular).  Another way to do this
could be to have a class (e.g. `Plant` friend a class `rPlant`) and
export just that.  Then the important guts of the model stay nice and
clean while the bits for R get kept off on the side.

Will require a lot of forwarding functions, but looking at the code,
that's mostly how it's done anyway.

But how deep do we do this for?  Looking at `Spline` it could easily
be done there, for sure.

However -- this causes its own issues.  At the moment, AdaptiveSpline
inherits from Spline (so that `init`, `eval`, `get_x`, `get_y`,
`get_xy` and `size` work on both classes).  We can't get these through
inheritance in R-wrapped versions because there would be no way of
initialising the base class sensibly (assuming that an `r_Spline`
would have a `Spline` object).

## Callbacks and functors

There are two basic ways that callbacks are being done --
AdaptiveSpline does it in an ad-hoc way, but the root finding (and the
quadrature approach) do the call backs with functors (I guess that the
ODE class is doing it a third way!).

Move everything that make sense over to functors.

Possible that one of `boost`s bind functions would help better here.

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

# Is the R code necessary?

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

# Efficiency issues

There will be many cases where the light environment is the same but
is recomputed (see notes above `Patch::ode_values_set`).  In addition,
there will be times where the light environment changes only because
of the death of a seedling (in which case only the points at the
seedling's height or lower need recomputing).  Ideally, it would be
nice to avoid as much of this calculation as possible.

However, the optimisation that seems immediately apparent to
`ode_values_set` will not work -- if the number of individuals has
changed this might report that it is unchanged even though it has.
However, it is possible that the ODE controller could know that we are
working with an OdeTarget and therefore nothing needs doing.

# Mutation and countable types

At the moment, everything is set up on the assumption of a small and
countable number of types.

If mutation is common and if every individual can have a different
strategy, the Species class disappears and we have a Patch of
individuals, each with their own strategy.

If mutation is rare enough, successful types may not have changed very
much and the Species idea might be useful enough to keep.

There may be a memory cost of having so many different individuals
with their own parameters.  I'm not sure that this will actually be a
problem in practice, but it could drastically increase the memory
requirements.  If it does, it's possible that there is a time/space
tradeoff where we store a set of "differences" from a base strategy,
where only a subset of traits actually change.  Then on lookup there
is a time cost for doing the two stage lookup -- possibly quite bad if
we have to go through a name lookup using strings or something.

# Unsorted issues

## Mortality

There are a couple of things weird about mortality:

* `mortality_rate` is the instantaneous death rate
* `mortality` is the integration of `mortality_rate` from some initial
  time to the present.
* `mortality_probability` is `1 - exp(-mortality_rate)` -- the
  probability that the plant has died between the initial time and the
  present (based on the 1 minus probability that exactly zero events
  happened with a mean number of events of `mortality_rate` from a
  Poisson distribution).

We only really need to use `survival_probability`, so I might remove
the `mortality_probability` case.

## Time

There is some consistency around the use of "time" in the model.

* The whole simulation has a single time axis.  This will start at
  zero and increase as the simulation runs.  Think of a calendar year,
  perhaps.  This is *the* time.
* A patch might need to know the time that it was born.  However,
  because I've not yet sorted out how disturbance will work in the
  stochastic version, this is not yet done.
* Individual `CohortTop` plants need to know when they are born *in
  the system time*.  This is done automatically by the
  `compute_initial_conditions`, which takes the system time from the
  environment.

## Consistency of get/set
All get/set pairs should be done with Rcpp's .property

`ode_values_set` -> `set_ode_values`
