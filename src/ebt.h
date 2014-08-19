// -*-c++-*-
#ifndef TREE_EBT_H_
#define TREE_EBT_H_

#include "patch.h"
#include "cohort_schedule.h"
#include "ode_solver.h"

namespace model {

class EBT : public ode::OdeTarget {
public:
  EBT(Parameters  p);
  EBT(Parameters *p);

  void run();
  std::vector<int> run_next();

  // * Output total seed rain calculation (not per capita)
  double seed_rain(size_t species_index) const;
  std::vector<double> seed_rains() const;
  std::vector<double> r_seed_rain_cohort(size_t species_index) const;
  std::vector<double> seed_rain_error(size_t species_index) const;
  std::vector<double> leaf_area_error(size_t species_index) const;

  double get_time() const;
  void reset();
  bool complete() const;

  // * ODE interface
  size_t ode_size() const;
  ode::iterator_const set_ode_values(double time_, ode::iterator_const it);
  ode::iterator       ode_values(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it)  const;

  // * R interface
  double r_seed_rain(size_t species_index) const;
  std::vector<double> r_seed_rain_error(size_t species_index) const;
  std::vector<double> r_leaf_area_error(size_t species_index) const;
  Patch<CohortTop> r_patch() const;
  CohortSchedule r_cohort_schedule() const;
  void r_set_cohort_schedule(CohortSchedule x);
  std::vector<double> r_ode_times() const;
  Parameters r_parameters() const;
  std::vector<double> r_times(size_t species_index) const;
  void r_set_times(std::vector<double> times, size_t species_index);
  void r_update_ode_times();

  Rcpp::List r_get_state() const;
  void r_set_state(Rcpp::List x);

private:
  std::vector<double> seed_rain_cohort(size_t species_index) const;

  Parameters::ptr parameters;
protected: // for now, at least (see EBTMutantRunner)
  Patch<CohortTop> patch;
private:
  ode::Solver<EBT> ode_solver;
  CohortSchedule schedule;
};

}

#endif
