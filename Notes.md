# Creating the shading environment

patch calls `compute_environment`, 
  - which calls `environment.compute_environment(f, height_max())`
    - where f = `compute_competiton` 
      -> compute_competition pass through species 
        - in species, compute_competition integrate density over nodes
          - todo: why not use QAG integrator here and anywhere else we do integration?
          -> in node, compute_competition = density * individual.compute_competition(height_); 
            -> individual passes straight through to strategy
              -> strategy calculates fraction of leaf above = k_I * area_leaf(height) * Q(z, height);
  - environment.compute_environment creates the spline of light
    - passes directly through to `canopy.compute_canopy(f_compute_competition, height_max)`
    - canopy_rescale_usually = false by default (check)

Can use the existing interpolator object, pass in vectors of x,y (of water)
see `AdaptiveInterpolator::update_spline()`

TODO:
- document classes QAG, tk etc
- rename `compute_competiton` to `precompute_direct_competiton`?
- rename `environment.compute_environment` ->  `environment.precompute_splines`?
  - similar: `canopy.compute_canopy` -> `shading.precompute_splines`?
- rename `canopy` -> `shading` or `direct_competiton_splines`? Canopy too specific?
- is environmet.rescale_environment used?
- `canopy::set_fixed_canopy` -> `canopy::set_fixed_values`

Environmentn.h already has vars.states (Internals object). Check Isaac not redefining these

# Using the shading environment

- In compute_rates, environment gets passed from patch -> Species -> Node -> Individual -> Strategy
  - FF16_strategy `net_mass_production_dt` calls `assimilation`
    - assimilation uses the object `quadrature::QAG integrator` to integrate the function `assimilation_leaf(environment.get_environment_at_height(z)) * q(z, height)` over z, where z is depth in the canopy. Based on the controls, the integration can be executed adaptively or not

TODO:

- How is Issac doing integration? What's happened to assimilation?


# Numerical methods

## Approximating functions with splines


## Numerical integration

### Gauss-Kronrod quadrature

The class `QAG` provides methods to use Gauss-Kronrod quadrature to numerically integration of a function, e.g. `double QAG::integrate(Function f, double a, double b)`. 

There are adaptive and non-adaptive methods. In the adaptive method, the difference between a Gauss quadrature rule and its Kronrod extension is used as an estimate of the approximation error, from which the need for more points is determined.

General background on the method is available at: https://en.wikipedia.org/wiki/Gaussian_quadrature
The QAG routine is adapted from the "QAG" algorithm in QUADPACK

The `QAG` class call class `QK`, which does the actual integration. This is code is ported from GSL.

- The integration has several "rules", defined in qk_rules.cpp. These include QK15, QK21, QK31, QK41, QK51. These allow for different numbers of points in the integration.





