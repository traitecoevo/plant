// -*-c++-*-
#ifndef NODE
#define NODE

#include <plant/models/environment.h>
#include <plant/numerics/gradient.h>
#include <plant/community/individual.h>
#include <odelia/ode_interface.hpp>
#include <plant/numerics/with_slope.h>
#include <limits> // std::numeric_limits

namespace plant {

template <typename T, typename E>
class Node {
public:
  using value_type = typename T::value_type;

  typedef T        strategy_type;
  typedef E        environment_type;
  typedef Individual<T,E> individual_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Node(strategy_type_ptr s);

  void compute_rates(const environment_type& environment, double pr_patch_survival);
  void compute_initial_conditions(const environment_type& environment, double pr_patch_survival, const value_type& birth_rate);

  // Wrapper to growth_rate_gradient for testing
  value_type r_growth_rate_gradient(const environment_type& environment);

  const value_type& height() const {return individual.state(HEIGHT_INDEX);}
  value_type compute_competition(const value_type& z) const;
  // The node's competition contribution and its vertical derivative, both
  // weighted by density. The first entry equals compute_competition(z) exactly.
  with_slope<value_type>
  compute_competition_and_slope(const value_type& z) const;
  value_type fecundity() const {
    return scaled_ ? offspring_produced_survival_weighted * density_scale
                   : offspring_produced_survival_weighted;
  }

  // Birth-date coordinate. The newest node and the boundary node are carried
  // per seed that establishes while the interval after the newest is still
  // open; their density and offspring are read through this factor, their
  // share of the establishment so far over the width it covers (see Species'
  // open_ode_size()). Species sets it on every state load; every other node
  // reads its state as is.
  void set_density_scale(const value_type& e) {
    scaled_ = true;
    density_scale = e;
  }
  void clear_density_scale() {
    scaled_ = false;
    density_scale = value_type(1.0);
  }
  value_type get_density_scale() const {
    return scaled_ ? density_scale : value_type(1.0);
  }
  // The next node has arrived, so this node's establishment is final: write it
  // into the state, where it stays. A node that nothing established on is
  // retired on the sentinels a cohort born with no density used to carry.
  void fold_establishment(const value_type& e) {
    scaled_ = false;
    density_scale = value_type(1.0);
    if (e > 0) {
      const value_type log_e = log(e);
      individual.set_state(MORTALITY_INDEX,
                           individual.state(MORTALITY_INDEX) - log_e);
      offspring_produced_survival_weighted *= e;
      set_log_density(log_density + log_e);
    } else {
      individual.set_state(MORTALITY_INDEX,
                           individual.state(MORTALITY_INDEX) +
                             establishment_failure_hazard);
      offspring_produced_survival_weighted = value_type(0.0);
      retire();
    }
  }

  // The survival factor offspring_produced_survival_weighted_dt multiplies. It
  // is zero rather than non-finite where the mortality integral has run away:
  // that happens only where the density is already too low to contribute, and a
  // non-finite factor here takes every rate reading this node with it.
  value_type survival_individual() const {
    const value_type s = exp(-individual.state(MORTALITY_INDEX));
    return util::is_finite(s) ? s : value_type(0.0);
  }

  // d(offspring_produced_survival_weighted_dt)/d(fecundity rate).
  value_type offspring_dt_dfecundity_rate(double pr_patch_survival) const {
    return survival_individual() * pr_patch_survival /
      pr_patch_survival_at_birth;
  }

  // The height growth rate this node was born at, i.e. |dh/dtau| exactly, at
  // birth. Current only for the boundary node (re-evaluated every step); frozen
  // at its own birth for an introduced one. Zero for a node loaded from an
  // exported state, which never ran compute_initial_conditions().
  double growth_rate_at_birth() const {return birth_growth_rate;}

  // Refresh only the birth date, leaving the rest of the bookkeeping alone.
  // Used for the not-yet-introduced boundary node, whose birth date is the
  // current time and so moves with every step; see
  // Species::set_new_node_birth_date().
  void set_introduction_time(double time) {node_introduction_time = time;}
  double introduction_time() const {return node_introduction_time;}
  double patch_density() const {return patch_density_at_birth;}
  double get_pr_patch_survival_at_birth() const {return pr_patch_survival_at_birth;}
  value_type get_log_density_rate() const {return log_density_dt;}

  // Restore birth bookkeeping for a node imported from an exported patch state,
  // without re-running compute_initial_conditions (which would overwrite the
  // loaded ODE state). pr_patch_survival_at_birth feeds the fecundity rate;
  // node_introduction_time and patch_density_at_birth feed lifetime-fitness
  // integrals. Required for a resumed run to reproduce the original trajectory.
  void set_birth_state(double time, double patch_density_in,
                       double pr_patch_survival) {
    node_introduction_time = time;
    patch_density_at_birth = patch_density_in;
    pr_patch_survival_at_birth = pr_patch_survival;
  }

  // Lifetime offspring of this node, weighted by the probability of
  // landing in a patch of the node's age and by survival during dispersal.
  value_type weighted_fecundity(const value_type& S_D) const {
    return fecundity() * patch_density_at_birth * S_D;
  }

  // A cohort carrying no density: absorbed into the top layer, or at exactly
  // zero density (born below the establishment threshold, or underflowed). Read
  // off the density itself, so a state loaded from a recording or an export
  // carries the answer with it and no record has to say which cohorts went.
  bool retired() const { return !(density > 0.0); }
  // Take the cohort's density to exactly zero, on the finite sentinel a cohort
  // that never established already carries: what absorption leaves behind, once
  // the density has moved into the top layer.
  void retire() {
    set_log_density(value_type(-establishment_failure_hazard));
  }

  // Unfortunate, but need a get_ here because of name shadowing...
  value_type get_log_density() const {return log_density;}
  // exp(log_density); can overflow to +Inf when the density equation runs
  // away (see Patch::check_finite_ode_state).
  // What the node contributes: its state, through the open interval's factor
  // where it has one.
  value_type get_density() const {
    return scaled_ ? density * density_scale : density;
  }
  // The same two values, read without recording: a check that only tests them
  // would otherwise put a statement on an active tape at every stage.
  double density_value() const {
    using odelia::util::to_passive;
    return scaled_ ? to_passive(density) * to_passive(density_scale)
                   : to_passive(density);
  }
  double log_density_value() const {
    return odelia::util::to_passive(log_density);
  }
  // log(get_density()), on the sentinel where that is zero.
  value_type effective_log_density() const {
    if (!scaled_) {
      return log_density;
    }
    return density_scale > 0 ? log_density + log(density_scale)
                             : value_type(-establishment_failure_hazard);
  }
  void set_log_density(const value_type& x) {
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
  template <typename It> It set_ode_state(It it);
  template <typename It> It ode_state(It it) const;
  template <typename It> It ode_rates(It it) const;
  template <typename It> It ode_aux(It it) const;
  template <typename It> It set_ode_aux(It it);

  static std::vector<std::string> ode_names() {
    std::vector<std::string> names = strategy_type::state_names();
    names.push_back("offspring_produced_survival_weighted");
    names.push_back("log_density");
    return names;
  }

  void resize_consumption_rates(int i) {
    individual.resize_consumption_rates(i);
  }

  value_type consumption_rate(int i) const {
    return individual.consumption_rate(i) * get_density();
  }

  individual_type individual;

  // Which coordinate this node's density is carried in. Resolves Control's
  // "auto" against the strategy's own requirement; see control.h.
  bool density_in_birth_date() const {
    return plant::density_in_birth_date(individual.control(),
                                        T::density_in_birth_date_default);
  }

  // The five this node carries beside the individual's. `density` is among them
  // although every load rewrites it: a write keeps whatever slot the value
  // already held.
  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, log_density, log_density_dt, density,
                              density_scale,
                              offspring_produced_survival_weighted,
                              offspring_produced_survival_weighted_dt,
                              individual);
  }

private:
  // The node's own two rates, from the individual's. Apart from compute_rates()
  // because compute_initial_conditions() has to seat the node's state between
  // the individual's rates and these, and BOTH OF THESE READ STATE THAT SEEDING
  // WRITES: the density rate through the guard below, the offspring rate through
  // survival_individual(). Taken before it, they describe the state the node
  // held a step ago.
  void compute_node_rates(const environment_type& environment,
                          double pr_patch_survival);

  // This is the gradient of growth rate with respect to height:
  value_type growth_rate_gradient(const environment_type& environment) const;

  value_type log_density;
  value_type log_density_dt;
  value_type density;  // exp(log_density), kept beside it for readers
  value_type density_scale{1.0};  // see set_density_scale()
  bool scaled_{false};
  value_type offspring_produced_survival_weighted;
  value_type offspring_produced_survival_weighted_dt;
  double pr_patch_survival_at_birth;

  // Recorded at introduction (see set_birth_state).
  double node_introduction_time;
  double patch_density_at_birth;
  // |dh/dtau| at birth; see growth_rate_at_birth().
  double birth_growth_rate;
};

template <typename T, typename E>
Node<T,E>::Node(strategy_type_ptr s)
  : individual(s),
    log_density(value_type(-std::numeric_limits<double>::infinity())),
    log_density_dt(0),
    density(0),
    offspring_produced_survival_weighted(0),
    offspring_produced_survival_weighted_dt(0),
    node_introduction_time(0),
    patch_density_at_birth(0),
    birth_growth_rate(0) {
}

template <typename T, typename E>
void Node<T,E>::compute_node_rates(const environment_type& environment,
                                   double pr_patch_survival) {
  // The coordinate branch lives in Individual::log_density_rate, which the step
  // recording reaches through this same call, so the recording and this path
  // cannot disagree about which coordinate they are on.
  //
  // ⚠️ A COHORT WITH NO DENSITY HAS NO LOG-DENSITY TRAJECTORY, and this guard is
  // what says so. On the HEIGHT coordinate `log_density_rate` is
  // `-d(growth)/d(height) - mortality`, which is non-zero for a dead cohort, so
  // without it the density DRIFTS OFF ITS FLOOR: exp(-establishment_failure_hazard)
  // is exactly zero, but integrating a positive rate up to -700 reaches 1e-305,
  // which is a cohort that does not exist acquiring a density. The -Inf this
  // sentinel replaced could not drift, so a finite one needs saying explicitly.
  //
  // ⚠️ A POSITIVE FINITE DENSITY, AND BOTH HALVES ARE LOAD-BEARING. `density > 0`
  // is what the floor needs once log(0)'s -Inf was replaced by a finite
  // sentinel; `is_finite(log_density)` is what an overflow to +Inf needs, which
  // passes `> 0` and is not a density anything can be carried along.
  log_density_dt = (util::is_finite(log_density) && density > 0.0)
                     ? individual.log_density_rate(environment)
                     : value_type(0.0);
  offspring_produced_survival_weighted_dt =
    individual.rate(FECUNDITY_INDEX) *
    offspring_dt_dfecundity_rate(pr_patch_survival);
}

template <typename T, typename E>
void Node<T,E>::compute_rates(const environment_type& environment,
                                double pr_patch_survival) {
  // A cohort with no density has no trajectory, and with retirement on it is
  // not integrated either: no physiology is evaluated for it, every rate it
  // reports is zero, and its individual stays where it froze. Gated on the top
  // layer (freezes_retired) rather than on the density alone, so a run that has
  // not opted in integrates its empty cohorts exactly as it always has.
  if (freezes_retired(individual.control()) && retired()) {
    individual.zero_rates();
    log_density_dt = value_type(0.0);
    offspring_produced_survival_weighted_dt = value_type(0.0);
    return;
  }
  individual.compute_rates(environment);
  compute_node_rates(environment, pr_patch_survival);
}

// On the height coordinate a new cohort's cumulative mortality starts at
// -log(establishment_probability), so the hazard of failing to establish is
// paid at birth, and its log_density at the log of the density it is born
// with. Where either is zero the state takes a finite sentinel instead:
// establishment_failure_hazard for the mortality, its negative for
// log_density. On the birth-date coordinate a cohort is born per seed that
// establishes, at hazard zero and the birth rate, and its species supplies the
// establishment over its cell (Species::set_open_scales()).
template <typename T, typename E>
void Node<T,E>::compute_initial_conditions(const environment_type& environment,
                                             double pr_patch_survival,
                                             const value_type& birth_rate) {
  pr_patch_survival_at_birth = pr_patch_survival;
  // Initialise strategy-specific initial states (e.g. TF24f's tracked position
  // at its optimum) before the first rates evaluation, so the birth growth rate uses
  // the initialised operating point rather than a default.
  individual.set_initial_states(environment);
  // On the birth-date coordinate the newborn's hazard is zero, and it is seated
  // before the rates so they read it, and not whatever the boundary node held
  // after the last evaluation -- which, after a retirement, is a hazard past
  // the gate, and a rate of zero.
  if (density_in_birth_date()) {
    individual.set_state(MORTALITY_INDEX, value_type(0.0));
  }
  // ⚠️ THE INDIVIDUAL'S RATES ONLY. The node's own two come after the initial
  // state below, because both read state it writes -- taken here, the density rate is
  // guarded against the density this node held BEFORE birth, which on one being
  // born for the first time is zero, and the rate it reports is that zero for
  // the whole of the node's first evaluation.
  individual.compute_rates(environment);

  // On the birth-date coordinate a node is born per seed that establishes, and
  // Species supplies the establishment through set_density_scale() while the
  // node's interval is open and folds it into the state when it closes.
  const value_type pr_estab = density_in_birth_date()
    ? value_type(1.0)
    : individual.establishment_probability_of_newborn(environment);
  // Held at establishment_failure_hazard on the zero arm rather than at the
  // -log(0) the equation gives; see the constant for what an infinite state
  // costs the reverse pass.
  individual.set_state("mortality",
                       pr_estab > 0 ? -log(pr_estab)
                                    : value_type(establishment_failure_hazard));
  // The birth-date axis of the node about to be introduced; Patch re-stamps
  // this with the exact introduction time as the node is pushed.
  node_introduction_time = environment.time;
  // dh/dtau = -g(H_0) at birth, so this is |dh/dtau| exactly, with no
  // differencing: characteristics are labelled by birth date, and a cohort born
  // an instant later starts an instant's growth behind. Recorded for every node
  // (a cheap read -- compute_rates() above has just set it, and the height
  // branch below needs it anyway), but only current for the *boundary* node,
  // which compute_initial_conditions() re-evaluates every step. For an
  // introduced node it is frozen at its own birth and so cannot serve as its
  // present-day Jacobian; see Species::height_jacobian().
  const value_type g = individual.rate(HEIGHT_INDEX);
  // A double diagnostic read back over the R boundary, so it takes the value.
  birth_growth_rate = odelia::util::to_passive(g);
  // The density itself, before its log: a recruit that cannot pay for itself
  // has an establishment probability of exactly zero on the arm below
  // threshold, so both this and its derivative are exactly zero there.
  const value_type density_at_birth =
    density_in_birth_date()
      ? birth_rate * pr_estab
      : (g > 0 ? birth_rate * pr_estab / g : value_type(0.0));
  // log() reaches -Inf from that zero through 0/0, so it records a non-finite
  // DERIVATIVE beside a correct value, and every rate reading the field this
  // node contributes to comes back non-finite. Only exp(log_density) is read
  // downstream, and it is continuous in both value and derivative here, so the
  // constructed constant is the same number carrying the zero the density has.
  //
  // ⚠️ AND THAT CONSTANT IS FINITE, for the reason the hazard's is: -Inf is a
  // STATE the patch carries for the rest of its life, and one no row can attach
  // to. `exp(-establishment_failure_hazard)` is exactly zero, so the density this
  // node contributes is the one log(0) gave it.
  set_log_density(density_at_birth > 0
                      ? log(density_at_birth)
                      : value_type(-establishment_failure_hazard));

  // A newborn per established seed whose hazard spends the establishment
  // failure hazard within the solver's shortest step survives no step the
  // solver can take: its survival is exactly zero in double precision by the
  // end of one. It is retired at birth, on the hazard a cohort that never
  // establishes carries, rather than integrated through a rate that steep.
  if (density_in_birth_date() &&
      !(individual.rate(MORTALITY_INDEX) * individual.control().ode_step_size_min <
        establishment_failure_hazard)) {
    individual.set_state(MORTALITY_INDEX,
                         value_type(establishment_failure_hazard));
    retire();
  }

  // Now that the density and the hazard are seated, the node's own rates.
  compute_node_rates(environment, pr_patch_survival);
}

template <typename T, typename E>
typename Node<T,E>::value_type
Node<T,E>::growth_rate_gradient(const environment_type& environment) const {
  return individual.growth_rate_gradient(environment);
}

// Wrapper to growth_rate_gradient for testing
template <typename T, typename E>
typename Node<T,E>::value_type
Node<T,E>::r_growth_rate_gradient(const environment_type& environment) {
  // We need to compute the physiological variables here, first, so
  // that reusing intervals works as expected.  This would ordinarily
  // be taken care of because of the calling order of
  // compute_rates / growth_rate_gradient.
  individual.compute_rates(environment);
  return growth_rate_gradient(environment);
}

template <typename T, typename E>
typename Node<T,E>::value_type
Node<T,E>::compute_competition(const value_type& height_) const {
  return get_density() * individual.compute_competition(height_);
}

template <typename T, typename E>
with_slope<typename Node<T,E>::value_type>
Node<T,E>::compute_competition_and_slope(const value_type& height_) const {
  const with_slope<value_type> fs =
    individual.compute_competition_and_slope(height_);
  const value_type d = get_density();
  return {d * fs.value, d * fs.slope};
}

// ODE interface -- note that the don't care about time in the node;
// only Patch and above does.
template <typename T, typename E>
template <typename It>
It Node<T,E>::set_ode_state(It it) {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    individual.set_state(i, *it++);
  }
  offspring_produced_survival_weighted = *it++;
  set_log_density(*it++);
  return it;
}
template <typename T, typename E>
template <typename It>
It Node<T,E>::ode_state(It it) const {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    util::write_iterator_scalar(it, individual.state(i));
  }
  util::write_iterator_scalar(it, offspring_produced_survival_weighted);
  util::write_iterator_scalar(it, log_density);
  return it;
}
template <typename T, typename E>
template <typename It>
It Node<T,E>::ode_rates(It it) const {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    util::write_iterator_scalar(it, individual.rate(i));
  }
  util::write_iterator_scalar(it, offspring_produced_survival_weighted_dt);
  util::write_iterator_scalar(it, log_density_dt);
  return it;
}

template <typename T, typename E>
template <typename It>
It Node<T,E>::ode_aux(It it) const {
  for (size_t i = 0; i < individual.aux_size(); i++) {
    util::write_iterator_scalar(it, individual.aux(i));
  }
  return it;
}

template <typename T, typename E>
template <typename It>
It Node<T,E>::set_ode_aux(It it) {
  return individual.set_ode_aux(it);
}


template <typename T, typename E>
Node<T,E> make_node(typename Node<T,E>::strategy_type s) {
  return Node<T,E>(make_strategy_ptr(s));
}

}

#endif /* NODE */
