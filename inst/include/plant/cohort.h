// -*-c++-*-
#ifndef PLANT_PLANT_COHORT_H_
#define PLANT_PLANT_COHORT_H_

#include <plant/environment.h>
#include <plant/gradient.h>
#include <plant/ode_interface.h>

namespace plant {

template <typename T>
class Cohort {
public:
  typedef T plant_type;
  typedef typename T::strategy_type   strategy_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Cohort(strategy_type_ptr s);

  void compute_vars_phys(const Environment& environment);
  void compute_initial_conditions(const Environment& environment);

  // * R interface (testing only, really)
  double r_growth_rate_gradient(const Environment& environment);

  double height() const {return plant.height();}
  double area_leaf_above(double z) const;
  double area_leaf() const;
  double fecundity() const {return seeds_survival_weighted;}

  // Unfortunate, but need a get_ here because of name shadowing...
  double get_log_density() const {return log_density;}
  void set_log_density(double x) {
    log_density = x;
    density = exp(log_density);
  }

  // ODE interface.
  //
  // NOTE: We are a time-independent model here so no need to pass
  // time in as an argument.  All the bits involving time are taken
  // care of by Environment for us.
  static size_t ode_size() {return 6;}
  ode::const_iterator set_ode_state(ode::const_iterator it);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;
  static std::vector<std::string> ode_names() {
    return std::vector<std::string>({"height", "mortality",
          "area_heartwood", "mass_heartwood",
          "seeds_survival_weighted", "log_density"});
  }

  plant_type plant;

private:
  // This is the gradient of growth rate with respect to height:
  double growth_rate_gradient(const Environment& environment) const;

  double log_density;
  double log_density_dt;
  double density; // hmm...
  double seeds_survival_weighted;
  double seeds_survival_weighted_dt;
  double pr_patch_survival_at_birth;
};

template <typename T>
Cohort<T>::Cohort(strategy_type_ptr s)
  : plant(s),
    log_density(R_NegInf),
    log_density_dt(0),
    density(0),
    seeds_survival_weighted(0),
    seeds_survival_weighted_dt(0),
    pr_patch_survival_at_birth(1) {
}

template <typename T>
void Cohort<T>::compute_vars_phys(const Environment& environment) {
  plant.compute_vars_phys(environment);

  // [eqn 22] Density per unit area of individuals with size 'h'.
  // EBT.md{eq:boundN}, see Numerical technique.
  // details.md, for details on translation from mass_leaf to height.
  // NOTE: This must be called *after* compute_vars_phys, but given we
  // need mortality_dt() that's always going to be the case.
  log_density_dt =
    - growth_rate_gradient(environment)
    - plant.mortality_dt();

  // EBT.md{eq:boundSurv}, see Numerical technique
  // survival_plant: converts from the mean of the poisson process (on
  // [0,Inf)) to a probability (on [0,1]).
  const double survival_patch = environment.patch_survival();
  double survival_plant = exp(-plant.mortality());
  if (!R_FINITE(survival_plant)) {
    // This is caused by NaN values in plant.mortality and log
    // density; this should only be an issue when density is so low
    // that we can throw these away.  I think that with smaller step
    // sizes this is better behaved too?
    survival_plant = 0.0;
  }

  seeds_survival_weighted_dt =
    plant.fecundity_dt() * survival_plant *
    survival_patch / pr_patch_survival_at_birth;
}

// NOTE: There will be a discussion of why the mortality rate initial
// condition is -log(germination_probability) in the documentation
// that Daniel is working out.
//
// NOTE: The initial condition for log_density is also a bit tricky, and
// defined on p 7 at the moment.
template <typename T>
void Cohort<T>::compute_initial_conditions(const Environment& environment) {

  compute_vars_phys(environment);
  // TODO: compute_vars_phys being computed again in germination_probability
  // removing here causing test failure, so figure out if can remove

  pr_patch_survival_at_birth = environment.patch_survival();
  const double pr_germ = plant.germination_probability(environment);
  // EBT.md{eq:boundSurv}
  plant.set_mortality(-log(pr_germ));
  // EBT.md{eq:boundN}
  const double g = plant.height_dt();
  const double seed_rain = environment.seed_rain_dt();
  // NOTE: log(0.0) -> -Inf, which should behave fine.
  set_log_density(g > 0 ? log(seed_rain * pr_germ / g) : log(0.0));

  // Need to check that the rates are valid after setting the
  // mortality value here (can go to -Inf and that requires squashing
  // the rate to zero).
  if (!R_FINITE(log_density)) {
    // Can do this at the same time that we do set_log_density, I think.
    log_density_dt = 0.0;
  }
  // NOTE: It's *possible* here that we need to set
  // plant.vars.mortality_dt to zero here, but I don't see that's
  // likely.
}

template <typename T>
double Cohort<T>::growth_rate_gradient(const Environment& environment) const {
  T p = plant;
  auto fun = [&] (double h) mutable -> double {
    return growth_rate_given_height(p, h, environment);
  };

  const Control& control = plant.control();
  const double eps = control.cohort_gradient_eps;
  if (control.cohort_gradient_richardson) {
    return util::gradient_richardson(fun, plant.height(), eps,
                                     control.cohort_gradient_richardson_depth);
  } else {
    return util::gradient_fd(fun, plant.height(), eps, plant.height_dt(),
                             control.cohort_gradient_direction);
  }
}

template <typename T>
double Cohort<T>::r_growth_rate_gradient(const Environment& environment) {
  // We need to compute the physiological variables here, first, so
  // that reusing intervals works as expected.  This would ordinarily
  // be taken care of because of the calling order of
  // compute_vars_phys / growth_rate_gradient.
  plant.compute_vars_phys(environment);
  return growth_rate_gradient(environment);
}

template <typename T>
double Cohort<T>::area_leaf_above(double height_) const {
  return density * plant.area_leaf_above(height_);
}

// TODO: Possibly push this logic into species and drop entirely from
// Cohort.
template <typename T>
double Cohort<T>::area_leaf() const {
  return area_leaf_above(0.0);
}

// ODE interface -- note that the don't care about time in the cohort;
// only Patch and above does.
template <typename T>
ode::const_iterator Cohort<T>::set_ode_state(ode::const_iterator it) {
  plant.set_height(*it++);
  plant.set_mortality(*it++);
  // skipping plant::fecunity, in lieu of seeds_survival_weighted
  plant.set_area_heartwood(*it++);
  plant.set_mass_heartwood(*it++);
  seeds_survival_weighted = *it++;
  set_log_density(*it++);
  return it;
}
template <typename T>
ode::iterator Cohort<T>::ode_state(ode::iterator it) const {
  *it++ = plant.height();
  *it++ = plant.mortality();
  *it++ = plant.area_heartwood();
  *it++ = plant.mass_heartwood();
  *it++ = seeds_survival_weighted;
  *it++ = log_density;
  return it;
}
template <typename T>
ode::iterator Cohort<T>::ode_rates(ode::iterator it) const {
  *it++ = plant.height_dt();
  *it++ = plant.mortality_dt();
  *it++ = plant.area_heartwood_dt();
  *it++ = plant.mass_heartwood_dt();
  *it++ = seeds_survival_weighted_dt;
  *it++ = log_density_dt;
  return it;
}

template <typename T>
Cohort<T> make_cohort(typename Cohort<T>::strategy_type s) {
  return Cohort<T>(make_strategy_ptr(s));
}

template <typename T>
double growth_rate_given_height(T& plant, double height,
                                const Environment& environment) {
  plant.set_height(height);
  plant.compute_vars_phys(environment, true);
  return plant.height_dt();
}

}

#endif
