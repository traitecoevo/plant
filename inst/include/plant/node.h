// -*-c++-*-
#ifndef NODE
#define NODE

#include <plant/environment.h>
#include <plant/gradient.h>
#include <plant/ode_interface.h>

namespace plant {

template <typename T, typename E>
class Node {
public:
  typedef T        strategy_type;
  typedef E        environment_type;
  typedef Individual<T,E> individual_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Node(strategy_type_ptr s);

  void compute_rates(const environment_type& environment, double pr_patch_survival);
  void compute_initial_conditions(const environment_type& environment, double pr_patch_survival, double birth_rate);

  // * R interface (testing only, really)
  double r_growth_rate_gradient(const environment_type& environment);

  double height() const {return individual.state(HEIGHT_INDEX);}
  double compute_competition(double z) const;
  double competition_effect() const;
  double fecundity() const {return offspring_produced_survival_weighted;}

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
  // +2 for log_density and offspring_production_dt
  static size_t ode_size() { return strategy_type::state_size() + 2; }
  size_t aux_size() const { return individual.aux_size(); }
  ode::const_iterator set_ode_state(ode::const_iterator it);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;
  ode::iterator       ode_aux(ode::iterator it) const;

  static std::vector<std::string> ode_names() {
    std::vector<std::string> names = strategy_type::state_names();
    names.push_back("offspring_produced_survival_weighted");
    names.push_back("log_density");
    return names;
  }

  void resize_consumption_rates(int i) {
    individual.resize_consumption_rates(i);
  }

  double consumption_rate(int i) const {
    return individual.consumption_rate(i) * density;
  }

  individual_type individual;

private:
  // This is the gradient of growth rate with respect to height:
  double growth_rate_gradient(const environment_type& environment) const;

  double log_density;
  double log_density_dt;
  double density; // hmm...
  double offspring_produced_survival_weighted;
  double offspring_produced_survival_weighted_dt;
  double pr_patch_survival_at_birth;
};

template <typename T, typename E>
Node<T,E>::Node(strategy_type_ptr s)
  : individual(s),
    log_density(R_NegInf),
    log_density_dt(0),
    density(0),
    offspring_produced_survival_weighted(0),
    offspring_produced_survival_weighted_dt(0) {
}

template <typename T, typename E>
void Node<T,E>::compute_rates(const environment_type& environment,
                                double pr_patch_survival) {
  individual.compute_rates(environment);

  // NOTE: This must be called *after* compute_rates, but given we
  // need mortality_dt() that's always going to be the case.
  log_density_dt =
    - growth_rate_gradient(environment)
    - individual.rate("mortality");

  // survival_individual: converts from the mean of the poisson process (on
  // [0,Inf)) to a probability (on [0,1]).
  double survival_individual = exp(-individual.state(MORTALITY_INDEX));
  if (!R_FINITE(survival_individual)) {
    // This is caused by NaN values in plant.mortality and log
    // density; this should only be an issue when density is so low
    // that we can throw these away.  I think that with smaller step
    // sizes this is better behaved too?
    survival_individual = 0.0;
  }

  offspring_produced_survival_weighted_dt =
    individual.rate("fecundity") * survival_individual *
    pr_patch_survival / pr_patch_survival_at_birth;
}

// NOTE: There will be a discussion of why the mortality rate initial
// condition is -log(establishment_probability) in the documentation
// that Daniel is working out.
//
// NOTE: The initial condition for log_density is also a bit tricky, and
// defined on p 7 at the moment.
template <typename T, typename E>
void Node<T,E>::compute_initial_conditions(const environment_type& environment,
                                             double pr_patch_survival, double birth_rate) {
  pr_patch_survival_at_birth = pr_patch_survival;
  compute_rates(environment, pr_patch_survival);

  const double pr_estab = individual.establishment_probability(environment);
  individual.set_state("mortality", -log(pr_estab));
  const double g = individual.rate("height");
  // NOTE: log(0.0) -> -Inf, which should behave fine.
  set_log_density(g > 0 ? log(birth_rate * pr_estab / g) : log(0.0));

  // Need to check that the rates are valid after setting the
  // mortality value here (can go to -Inf and that requires squashing
  // the rate to zero).
  if (!R_FINITE(log_density)) {
    // Can do this at the same time that we do set_log_density, I think.
    log_density_dt = 0.0;
  }
  // NOTE: It's *possible* here that we need to set
  // individual.vars.mortality_dt to zero here, but I don't see that's
  // likely.
}

template <typename T, typename E>
double Node<T,E>::growth_rate_gradient(const environment_type& environment) const {
  individual_type p = individual;
  auto fun = [&] (double h) mutable -> double {
    return p.growth_rate_given_height(h, environment);
  };

  // std::cout<<"\t  grad\t";
  const Control& control = individual.control();
  const double eps = control.node_gradient_eps;
  if (control.node_gradient_richardson) {
    return util::gradient_richardson(fun,  individual.state(HEIGHT_INDEX), eps,
                                     control.node_gradient_richardson_depth);
  } else {
    return util::gradient_fd(fun, individual.state(HEIGHT_INDEX), eps, individual.rate("height"),
                             control.node_gradient_direction);
  }
}

template <typename T, typename E>
double Node<T,E>::r_growth_rate_gradient(const environment_type& environment) {
  // We need to compute the physiological variables here, first, so
  // that reusing intervals works as expected.  This would ordinarily
  // be taken care of because of the calling order of
  // compute_rates / growth_rate_gradient.
  individual.compute_rates(environment);
  return growth_rate_gradient(environment);
}

template <typename T, typename E>
double Node<T,E>::compute_competition(double height_) const {
  return density * individual.compute_competition(height_);
}

template <typename T, typename E>
double Node<T,E>::competition_effect() const {
  return compute_competition(0.0);
}

// ODE interface -- note that the don't care about time in the node;
// only Patch and above does.
template <typename T, typename E>
ode::const_iterator Node<T,E>::set_ode_state(ode::const_iterator it) {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    individual.set_state(i, *it++);
  }
  offspring_produced_survival_weighted = *it++;
  set_log_density(*it++);
  return it;
}
template <typename T, typename E>
ode::iterator Node<T,E>::ode_state(ode::iterator it) const {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    *it++ = individual.state(i);
  }
  *it++ = offspring_produced_survival_weighted;
  *it++ = log_density;
  return it;
}
template <typename T, typename E>
ode::iterator Node<T,E>::ode_rates(ode::iterator it) const {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    *it++ = individual.rate(i);
  }
  *it++ = offspring_produced_survival_weighted_dt;
  *it++ = log_density_dt;
  return it;
}

template <typename T, typename E>
ode::iterator Node<T,E>::ode_aux(ode::iterator it) const {
  for (size_t i = 0; i < individual.aux_size(); i++) {
    *it++ = individual.aux(i);
  }
  return it;
}

template <typename T, typename E>
Node<T,E> make_node(typename Node<T,E>::strategy_type s) {
  return Node<T,E>(make_strategy_ptr(s));
}

}

#endif /* NODE */
