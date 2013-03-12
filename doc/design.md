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

I've been using Java style file names - switch to lowercase names with
camel->underscore.

Deciding if to go with class variables with underscore (google style
guide).

# Consistency

init vs contructors
reset, clear, etc.

