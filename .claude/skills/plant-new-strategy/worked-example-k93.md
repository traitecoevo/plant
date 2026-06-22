# Worked example: implementing Kohyama 1993 (K93) as a new strategy

This is the reference walkthrough for the **biology** step of `plant-new-strategy`
(originally `vignettes/strategy_new.Rmd`). It shows how the shipped K93 model was
built by scaffolding from FF16 and then filling in a simple size-structured
model. K93 already exists in the package — read its files
(`inst/include/plant/models/k93_strategy.h`, `src/k93_strategy.cpp`, `R/k93.R`)
alongside this guide.

K93 tracks growth, reproduction and mortality as functions of stem size and the
cumulative basal area of taller competitors. Disturbance (not Kohyama's original
formulation) supplies background mortality; dispersal/establishment are set to 1.

## Model equations

- **Competition (basal area above size x):** `B(x) = (π/4) ∫ y² f(y) dy` for
  `y ≥ x` — the basal area of all individuals taller than `x`.
- **Growth (eqn 10):** `G(x) = x · (b₀ − b₁·ln x − b₂·B(x))`, floored at 0.
- **Fecundity (eqn 12):** `R(x) = d₀ · basal_area(x) · exp(−d₁·B(x))`.
- **Mortality (eqn 11):** `μ(x) = −c₀ + c₁·B(x)`, floored at 0.

## Strategy source (`src/k93_strategy.cpp`)

```cpp
// [eqn 10] Growth
double K93_Strategy::size_dt(double size, double cumulative_basal_area) const {
  double growth = size * (b_0 - b_1 * log(size) - b_2 * cumulative_basal_area);
  return (growth < 0.0) ? 0.0 : growth;
}

// [eqn 12] Reproduction
double K93_Strategy::fecundity_dt(double size, double cumulative_basal_area) const {
  double basal_area = size_to_basal_area(size);
  return d_0 * basal_area * exp(-d_1 * cumulative_basal_area);
}

// [eqn 11] Mortality
double K93_Strategy::mortality_dt(double cumulative_basal_area,
                                  double cumulative_mortality) const {
  // If survival probability has already collapsed (latency = Inf) the rate
  // calc breaks; returning 0 gives the correct behaviour.
  if (R_FINITE(cumulative_mortality)) {
    double mu = -c_0 + c_1 * cumulative_basal_area;
    return (mu > 0) ? mu : 0.0;
  }
  return 0.0;
}

// map size -> basal area, and the competition contribution
double K93_Strategy::size_to_basal_area(double size) const {
  return M_PI / 4 * pow(size, 2);
}

double K93_Strategy::compute_competition(double z, double size) const {
  // competition felt only by individuals shorter than the focal size
  return size_to_basal_area(size) * Q(z, size);
}
```

## Environment mapping (`*_environment.h`)

The environment turns the cumulative competition into a light-availability
spline in `[0,1]` (queried per individual as an O(1) lookup):

```cpp
void compute_environment(Function f_compute_competition, double height_max, bool rescale) {
  const double lower_bound = 0.0;
  double upper_bound = height_max;
  auto f_light_availability = [&](double height) -> double {
    return exp(-k_I * f_compute_competition(height));
  };
  environment_interpolator =
    environment_generator.construct(f_light_availability, lower_bound, upper_bound);
}
```

## Wiring rates (`compute_rates`)

Back-transform the `[0,1]` light value to cumulative basal area, then set the
ODE rates by integer state index:

```cpp
void K93_Strategy::compute_rates(const K93_Environment& environment, Internals& vars) {
  double height = vars.state(HEIGHT_INDEX);

  double competition = environment.get_environment_at_height(height);
  double cumulative_basal_area = -log(competition) / environment.k_I;
  if (!util::is_finite(cumulative_basal_area)) {
    util::stop("Environmental interpolation out of bounds");
  }

  vars.set_rate(HEIGHT_INDEX,    size_dt(height, cumulative_basal_area));
  vars.set_rate(FECUNDITY_INDEX, fecundity_dt(height, cumulative_basal_area));
  vars.set_rate(MORTALITY_INDEX, mortality_dt(cumulative_basal_area,
                                              vars.state(MORTALITY_INDEX)));
}
```

## Header declarations (`*_strategy.h`)

The subclass declares its parameters, the special methods the core requires, and
nuisance fixtures:

```cpp
class K93_Strategy : public Strategy<K93_Environment> {
public:
  typedef std::shared_ptr<K93_Strategy> ptr;
  K93_Strategy();

  // --- the two methods the core requires ---
  void   compute_rates(const K93_Environment& environment, Internals& vars);
  double compute_competition(double z, double size) const;

  // --- parameters ---
  double height_0;             // initial seedling size (dbh cm)
  double b_0, b_1, b_2;        // growth: intercept, asymptote, suppression
  double d_0, d_1;             // reproduction: rate, suppression
  double c_0, c_1;             // mortality: intercept, suppression

  // --- methods ---
  double size_to_basal_area(double size) const;
  double size_dt(double size, double cumulative_basal_area) const;
  double fecundity_dt(double size, double cumulative_basal_area) const;
  double mortality_dt(double cumulative_basal_area, double cumulative_mortality) const;

  // --- fixtures ---
  double S_D = 1.0;            // dispersal survival (required by scm.h)
  double eta = 12;             // competition smoothing

  double establishment_probability(const K93_Environment& environment) { return 1.0; }

  // smoothing of the competition edge
  double Q(double z, double size) const {
    if (z > size) return 0.0;
    const double tmp = 1.0 - pow(z / size, eta);
    return tmp * tmp;
  }

  void prepare_strategy();     // set constants
  std::string name;
};
```

Keep `state_names()`/`state_size()` and the `*_INDEX` `constexpr`s in sync (see
agents.md §11/§12). K93 uses `height, mortality, fecundity`.

## RcppR6 yml block

The scaffolder copies the template's `<Name>_Strategy:` block — trim the `list:`
to the parameters the model actually exposes:

```yaml
K93_Strategy:
  name_cpp: "plant::K93_Strategy"
  roxygen: |
    Strategy parameters that tune various aspects of the biological model.
    @title Strategy parameters
    @param ...,values Values to initialise the struct with (either as
    variadic arguments, or as a list, but not both).
    @export
  list:
    - height_0: double
    - b_0: double
    - b_1: double
    - b_2: double
    - c_0: double
    - c_1: double
    - d_0: double
    - d_1: double
    - control: "plant::Control"
```

(`K93_Environment` needed no changes beyond the scaffold.)

## R hyperparameters (`R/k93.R`)

`make_<Name>_hyperpar` returns a closure mapping a trait matrix → strategy
parameters. K93's defaults (Kohyama 1993, Table/§ values):

```r
make_K93_hyperpar <- function(
    b_0 = 0.059,    # growth intercept (yr⁻¹)
    b_1 = 0.012,    # growth asymptote (yr⁻¹·(ln cm)⁻¹)
    b_2 = 0.00041,  # growth suppression (m²·cm⁻²·yr⁻¹)
    c_0 = 0.008,    # mortality intercept (yr⁻¹)
    c_1 = 0.00044,  # mortality suppression (m²·cm⁻²·yr⁻¹)
    d_0 = 0.00073,  # recruitment rate (cm²·yr⁻¹)
    d_1 = 0.044) {  # recruitment suppression (m²·cm⁻²)
  # assert each is scalar ...
  function(m, s, filter = TRUE) m   # K93 sets params directly; no trade-offs
}
```

## Build and run

```r
# make rebuild   # at the shell, after editing the interface
p0 <- scm_base_parameters("K93")
p1 <- expand_parameters(trait_matrix(0.0825, "b_0"), p = p0)
p1$birth_rate <- 20
res <- run_scm(p1, collect = TRUE, refine_schedule = TRUE)
# res$species is a tidy tibble: time, node, height, mortality, fecundity, ...
```

K93 was tuned with a ~200 yr disturbance interval (vs FF16's ~30 yr).
