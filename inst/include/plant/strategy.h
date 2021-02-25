// -*-c++-*-
#ifndef PLANT_PLANT_STRATEGY_H_
#define PLANT_PLANT_STRATEGY_H_

#include <memory>
#include <plant/control.h>
#include <plant/internals.h>
#include <plant/uniroot.h>
#include <RcppCommon.h> // NA_REAL
#include <plant/uniroot.h>


namespace plant {

template <typename E> 
class Strategy {
public:
  typedef E             environment_type;
  typedef std::shared_ptr<Strategy> ptr;

  // update this when the length of state_names changes
  static size_t state_size ();
  // update this when the length of aux_names changes
  size_t aux_size ();

  static std::vector<std::string> state_names();

  std::vector<std::string> aux_names();

  // TODO : expose this so can access state_names directly
  // In previous attempt couldn't get it to run
  // static std::vector<std::string> state_names() { return strategy_type::state_names(); }
  // the index of variables in the internals extra vector
  std::map<std::string, int> state_index; 
  std::map<std::string, int> aux_index; 


  bool collect_all_auxillary;

  void refresh_indices();

  double competition_effect(double size) const;

  double competition_effect_state(Internals& vars);

  void compute_rates(const environment_type& environment, bool reuse_intervals,
                Internals& vars);

  void update_dependent_aux(const int index, Internals& vars);

  double net_mass_production_dt(const environment_type& environment,
                                double size, double competition_effect_,
                                bool reuse_intervals=false);

  double establishment_probability(const environment_type& environment);

  double fecundity_dt(double net_mass_production_dt,
                      double fraction_allocation_reproduction) const;

  double mortality_dt(double productivity_area, double cumulative_mortality) const;

  double compute_competition(double z, double size) const;

  double initial_size(void) const;

  double size_0;

  // Every Strategy needs a set of Control objects -- these govern
  // things to do with how numerical calculations are performed,
  // rather than the biological control that this class has.
  Control control;

  std::string name;
};


}

#endif
