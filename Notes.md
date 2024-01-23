# Creating the shading environment

patch calls `compute_environment`, 
  - which calls `environment.compute_environment(f, height_max())`
    - where f = `compute_competiton` 
      -> compute_competition pass through species 
        - in species, compute_competition integrate density over nodes
          - todo: why not use QAG integrator here?
          -> in node, compute_competition = density * individual.compute_competition(height_); 
            -> individual passes straight through to strategy
              -> strategy calculates fraction of leaf above = k_I * area_leaf(height) * Q(z, height);
  - environment.compute_environment creates the spline of light
    - passes directly through to `canopy.compute_canopy(f_compute_competition, height_max)`
    - canopy_rescale_usually = false by default

Can use the existing interpolator object, pass in vectors of x,y (of water)
see `AdaptiveInterpolator::update_spline()`

TODO:
- rename `compute_competiton` to `compute_direct_competiton`?
- rename `environment.compute_environment` ->  `environment.precompute_splines`?
  - similar: `canopy.compute_canopy` -> `shading.precompute_splines`?
- rename `canopy` -> `shading` or `direct_competiton_splines`? Canopy too specific?
- is environmet.rescale_environment used?
- `canopy::set_fixed_canopy` -> `canopy::set_fixed_values`

Environmentn.h already has vars.states (Internals object). Check Isaac not redefining these

# Using the shading environment

- In compute_rates, environment gets passed from patch -> Species -> Node -> Individual -> Strategy
  - FF16_strategy `net_mass_production_dt` calls `assimilator.assimilate(environment, height,
                                            area_leaf_, reuse_intervals)`
    - assimilator has object `quadrature::QAG integrator`, which does the integrations, adaptively or not;

TODO:

- would it be more useful to have assimilator take function for leaf ptss and weights, and average it, i.e. store leaf model, canopy in FF16? 
- How is Issac doing integration? What's happened to assimilator?
