# Notes on the EBT version of the model

The method `EBT_Base::calculate_rates_of_change` seems to be the
equivalent of the derivatives function, or at least one of them; it's
called within `EBT_Base::take_step`.  The other half looks to be
`EBT_Base::update_states`; the step looks like this

```
calculate_rates_of_change();
for each part of the RK tableau:
  compute light environment with make_env_spline()
  calculate_rates_of_change();
  update_states();
```

## State information:

* `mu`: mass of average individual
* `lam`: total number of individuals
* `coh_R`: seed output
* `bound`: upper bound (in leaf size (?))
* `bound_n0`: density of largest individuals at first time
* `bound_nt`: density at current time
* `bound_S`: survival at cohort boundary
* `bound_R`: cumulative (lifetime) seed production at cohort boundary
* `bottom_coh_mass`: mass of the bottom cohort
* `height`: mean height, given mean mass

## `calculate_rates_of_change`:

First, compute disturbance rate at current time (= patch age).  The
rest of the method is a loop over all species.  For each species:

1. Zero key variables (`Seed_Rain_out`, `NPP`, `GPP`, `Turnover`,
   `Resp`).  (*I don't think we track all of these, but could probably
   arrange to get them if we really cared*)
2. Looping over the cohorts from smallest to largest:
   a. Let `size` be the size of the upper bound
   b. Compute all size and physiological variables for something this
	  (using `update_states`)
   c. Collect into `d_bound`, `d_bound_S`, and `d_bound_R` the growth
      rate, *negative* mortality rate and fecundity rate, with
      `d_bound_R` computed as `f * exp(bound_S) * w(birth, t)`
   d. Compute `d_bound_dt` by tweaking the height slightly (a finite
      difference approximation of fixed size).  If `g'` is the growth
      rate for a cohort of size `size - ds`, then
	  `d_bound_dt = -(d + (g - g')/ds)` (where `d` is the mortality
      rate).
   e. Set `d_lambda` equal to `-d * lambda` (decrease by the
      per-capita mortality rate multiplied by the population size).
   f. If this is the bottom cohort, the rate in change in the bottom
      mass is `d_bottom_coh_mass` is `g*lam - d * bottom_coh_mas`
	  (*I would have thought that it would stay at the seed leaf mass
      size?*)
   g. The fecundity rate `coh_R` is `f * lam * patch_age_density_freq(time)`
   h. Then species-level rates are incremented: 
      - `Seed_rain_out` by `d_r`
	  - `GPP` by `A * lam`
	  - `NPP` by `N * lam`
	  - `Turnover` by `T * lam`
	  - `Resp` by `R * lam`

During this, there is special treatment for the bottom (smallest)
cohort.  Before doing anything, we set `S_germ` equal to the
germination probability in the current light environment, and set
`bound_S` to the *log* of this (*why the log?*).

## Other questions

Discuss the "seed rain" more clearly - what is it doing?

## Level of control

The "patch weight" and "time" calculation comes in a couple of places;

* Calculating the disturbance rate (possibly never used)
* Computing `patch_weight` for `d_bound_R`
* Computing `patch_age_density_freq` for `d_R`

Probably the simplest way to deal with these would be to set them at
the same time as we compute the light environment; in some ways this
is really part of the environment.

The other way of dealing with these is to do nothing with them until
later on (at the species or patch level).  They apparently only feed
in to the calculations around fecundity, so that might not be too
hard.

# Notes

We compute the cohort weight as `cohort_weight(time_at_birth,
current_time`; this could be a method of the cohort, so that the a
method `weight(double current_time)` returns the weight of that
cohort.  However, that would require that a cohort knew about the
disturbance regime (wrong level).  A patch could concievably know when
a cohort was born.
