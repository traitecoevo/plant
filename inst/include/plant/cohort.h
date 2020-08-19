// -*-c++-*-
#ifndef PLANT_PLANT_COHORT_H_
#define PLANT_PLANT_COHORT_H_

#include <plant/environment.h>
#include <plant/gradient.h>
#include <plant/ode_interface.h>

namespace plant {

template <typename T, typename E>
class Cohort {
public:
  typedef T        strategy_type;
  typedef E        environment_type;
  typedef Individual<T,E> individual_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Cohort(strategy_type_ptr s);

  void compute_rates(const environment_type& environment);
  void compute_initial_conditions(const environment_type& environment);

  // * R interface (testing only, really)
  double r_growth_rate_gradient(const environment_type& environment);

  double height() const {return plant.state(HEIGHT_INDEX);}
  double compute_competition(double z) const;
  double competition_effect() const;
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
  static size_t ode_size() { return strategy_type::state_size() + 2; }
  ode::const_iterator set_ode_state(ode::const_iterator it);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;

  static std::vector<std::string> ode_names() {
    std::vector<std::string> plant_names = strategy_type::state_names();
    plant_names.push_back("seeds_survival_weighted");
    plant_names.push_back("log_density");
    return plant_names;
  }

  individual_type plant;

private:
  // This is the gradient of growth rate with respect to height:
  double growth_rate_gradient(const environment_type& environment) const;

  double log_density;
  double log_density_dt;
  double density; // hmm...
  double seeds_survival_weighted;
  double seeds_survival_weighted_dt;
  double pr_patch_survival_at_birth;
};

template <typename T, typename E>
Cohort<T,E>::Cohort(strategy_type_ptr s)
  : plant(s),
    log_density(R_NegInf),
    log_density_dt(0),
    density(0),
    seeds_survival_weighted(0),
    seeds_survival_weighted_dt(0),
    pr_patch_survival_at_birth(1) {
}

template <typename T, typename E>
void Cohort<T,E>::compute_rates(const environment_type& environment) {
  plant.compute_rates(environment);

  // NOTE: This must be called *after* compute_rates, but given we
  // need mortality_dt() that's always going to be the case.
  log_density_dt =
    - growth_rate_gradient(environment)
    - plant.rate("mortality");

  // survival_plant: converts from the mean of the poisson process (on
  // [0,Inf)) to a probability (on [0,1]).
  const double survival_patch = environment.patch_survival();
  double survival_plant = exp(-plant.state(MORTALITY_INDEX));
  if (!R_FINITE(survival_plant)) {
    // This is caused by NaN values in plant.mortality and log
    // density; this should only be an issue when density is so low
    // that we can throw these away.  I think that with smaller step
    // sizes this is better behaved too?
    survival_plant = 0.0;
  }

  seeds_survival_weighted_dt =
    plant.rate("fecundity") * survival_plant *
    survival_patch / pr_patch_survival_at_birth;
}

// NOTE: There will be a discussion of why the mortality rate initial
// condition is -log(establishment_probability) in the documentation
// that Daniel is working out.
//
// NOTE: The initial condition for log_density is also a bit tricky, and
// defined on p 7 at the moment.
template <typename T, typename E>
void Cohort<T,E>::compute_initial_conditions(const environment_type& environment) {
  compute_rates(environment);

  pr_patch_survival_at_birth = environment.patch_survival();
  const double pr_germ = plant.establishment_probability(environment);
  plant.set_state("mortality", -log(pr_germ));
  const double g = plant.rate("height");
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

template <typename T, typename E>
double Cohort<T,E>::growth_rate_gradient(const environment_type& environment) const {
  individual_type p = plant;
  auto fun = [&] (double h) mutable -> double {
    return p.growth_rate_given_height(h, environment);
  };

  const Control& control = plant.control();
  const double eps = control.cohort_gradient_eps;
  if (control.cohort_gradient_richardson) {
    return util::gradient_richardson(fun,  plant.state(HEIGHT_INDEX), eps,
                                     control.cohort_gradient_richardson_depth);
  } else {
    return util::gradient_fd(fun, plant.state(HEIGHT_INDEX), eps, plant.rate("height"),
                             control.cohort_gradient_direction);
  }
}

template <typename T, typename E>
double Cohort<T,E>::r_growth_rate_gradient(const environment_type& environment) {
  // We need to compute the physiological variables here, first, so
  // that reusing intervals works as expected.  This would ordinarily
  // be taken care of because of the calling order of
  // compute_rates / growth_rate_gradient.
  plant.compute_rates(environment);
  return growth_rate_gradient(environment);
}

template <typename T, typename E>
double Cohort<T,E>::compute_competition(double height_) const {
  return density * plant.compute_competition(height_);
}

template <typename T, typename E>
double Cohort<T,E>::competition_effect() const {
  return compute_competition(0.0);
}

// ODE interface -- note that the don't care about time in the cohort;
// only Patch and above does.
template <typename T, typename E>
ode::const_iterator Cohort<T,E>::set_ode_state(ode::const_iterator it) {
  for (int i = 0; i < plant.ode_size(); i++) {
    plant.set_state(i, *it++);
  }
  seeds_survival_weighted = *it++;
  set_log_density(*it++);
  return it;
}
template <typename T, typename E>
ode::iterator Cohort<T,E>::ode_state(ode::iterator it) const {
  for (int i = 0; i < plant.ode_size(); i++) {
    *it++ = plant.state(i);
  }
  *it++ = seeds_survival_weighted;
  *it++ = log_density;
  return it;
}
template <typename T, typename E>
ode::iterator Cohort<T,E>::ode_rates(ode::iterator it) const {
  for (int i = 0; i < plant.ode_size(); i++) {
    *it++ = plant.rate(i);
  }
  *it++ = seeds_survival_weighted_dt;
  *it++ = log_density_dt;
  return it;
}

template <typename T, typename E>
Cohort<T,E> make_cohort(typename Cohort<T,E>::strategy_type s) {
  return Cohort<T,E>(make_strategy_ptr(s));
}

}

#endif
