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
  void run_next();

  // * Fitness calculation
  double fitness(size_t species_index) const;
  std::vector<double> fitnesses() const;
  std::vector<double> r_fitness_cohort(size_t species_index) const;
  std::vector<double> fitness_error(size_t species_index) const;
  std::vector<double> leaf_area_error(size_t species_index) const;

  double get_time() const;
  void reset();

  // * ODE interface
  size_t ode_size() const;
  ode::iterator_const set_ode_values(double time_, ode::iterator_const it);
  ode::iterator       ode_values(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it)  const;

  // * R interface
  double r_fitness(size_t species_index) const;
  std::vector<double> r_fitness_error(size_t species_index) const;
  std::vector<double> r_leaf_area_error(size_t species_index) const;
  Patch<CohortTop> r_patch() const;
  CohortSchedule r_cohort_schedule() const;
  void r_set_cohort_schedule(CohortSchedule x);
  std::vector<double> r_ode_times() const;
  Parameters r_parameters() const;
  std::vector<double> r_times(size_t species_index) const;
  void r_set_times(std::vector<double> times, size_t species_index);

  Rcpp::List r_get_state() const;
  void r_set_state(Rcpp::List x);

private:
  std::vector<double> fitness_cohort(size_t species_index) const;

  Parameters::ptr parameters;
  Patch<CohortTop> patch;
  ode::Solver<Patch <CohortTop> > ode_solver;
  CohortSchedule schedule;
};

}

#endif
