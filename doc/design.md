# Overall

This represents our "version 1" of the model.  In it we will have a
deterministic model and a stochastic individual based simulation.
However, the scope is very much the same as the original "EBT" code.

 - no explit space at the within patch level (but possibly at the
   between patch level).
 - one size dimension (but possibly more for the individual based
   model).
 - large number of possibly varying traits
 
# Shared components

## Individuals

The "individuals" in the individual based model are very similar to
cohorts in the deterministic model; they have a state (mass, leaf
area, etc) and are capable of computing their physiological rates
given a light environment.

However, a "cohort" will actually contain two individuals (top and
middle), and as a result they will provide different numbers of rates
back to the ODE solvers, etc.

## Species 

Above the level of species, there are two places that "species"
appears as a concept:

1. as a structure/class of parameters, corresponding to a different
   strategy.
2. as a book-keeping device, collecting individuals within a
   population
   
## Patch

A patch contains multiple species.  In the deterministic model, there
will be just one patch.

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

So, define a class member `standalone`, which will be false normally
but true if a `Plant` is created without a `Strategy` pointer.  In
that case, it will either create the default `Strategy` or take a copy
of a provided `Strategy` object, and be responsible for cleaning up on
deletion.  Copies involving a standalone `Plant` will *copy* the
`Parameter` object and create a new "standalone" `Plant`.

## Constants

There are probably still some constants that are settable at the level
of `Strategy`.

For example, `(c_acc * s)` appear together and only together in
computing `fecundity_rate`.  As such they are entirely colinear.

In theory `Pi_0` could be added here too, as the mortality during
dispersal (as was done in the EBT) but that does mix ideas a bit too
much.

In computing `mortality`, the value
`strategy->c_d0 * exp(-strategy->c_d1*strategy->rho)`
is constant.

There are others in `compute_vars_size`:
  - `rho / theta * a1 * eta_c` (`sapwood`)
  - `rho         * a2 * eta_c` (`heartwood`)

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

## ODE solver

We need an ODE solver that can cope with changing dimensions, as both
the individual-based and deterministic solvers have rapidly changing
dimensions (and the GSL solvers will have to reallocate every step,
which will become annoying).

## Splines

We use splines for the light environment, but I have some code I'm
quite happy with that works here.

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


# Consistency

init vs contructors
reset, clear, etc.

